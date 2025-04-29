classdef ncep < handle
    properties (GetAccess = public, SetAccess = protected) 
        model (1,1) string;
        product (1,1) string;
        cycle (1,1) datetime;
        depth (1,1) duration;
        local (1,1) string = missing;
        origin (1,1) string;
        inventory table;
    end

    methods
        function ref = ncep(model, product, cycle, depth)
            arguments
                model (1,1) string;
                product (1,1) string;
                cycle (1,1) datetime {mustBeFinite, ncep.shouldHaveTimezone};
                depth (1,1) duration {mustBeFinite, mustBeNonempty};
            end
            
            if isempty(cycle.TimeZone)
                warning("ncep:timezone", ...
                    "Model cycle does not have assigned time zone. Defaulting to UTC.");
            end
            cycle.TimeZone = "UTC";
            cycle.Format = "dd-MMM-yyyy HH:mm z";
            cycle_day = dateshift(cycle, "start", "day");
            cycle_time = cycle - cycle_day;

            [~, source] = ncep.list(model, product, cycle_time, depth);
            ref.model = model;
            ref.product = product;
            ref.cycle = cycle;
            ref.depth = depth;

            ref.origin = ncep.sprintnf(source, ...
                product = ref.product, cycle = ref.cycle.Hour, fcst = hours(ref.depth), ...
                year = ref.cycle.Year, month = ref.cycle.Month, day = ref.cycle.Day);
            ref.inventory = ncep.read_ncep_index(ref.origin + ".idx");
            ref.inventory = ncep.patch_index(ref.inventory);
        end

        function [data, crs] = read(refs, params)
            arguments
                refs (1, :) ncep;
                params.lats (1,2) double = [-Inf Inf];
                params.lons (1,2) double = [-Inf Inf];
                params.field (1, :) pattern = wildcardPattern;
                params.layer (1, :) pattern = wildcardPattern;
            end

            [params, lats] = ncep.strip(params, "lats");
            [params, lons] = ncep.strip(params, "lons");
            for i_ref = 1:length(refs)
                refs(i_ref).download(params);
            end

            % create template from first object (though they should all be the same)
            template = refs(1).inventory;
            params.band = @isfinite; % filter to downloaded data (unset bands are NaN)
            messages = ncep.search(template, params);
            template = template(messages, :) ;

            output_fields = unique(template.field);
            output_layers = unique(template.layer);
            times = [refs.cycle] + [refs.depth];

            info = georasterinfo(refs(1).local);
            [~, data_indices, axes] = ncep.crop_raster(info.RasterReference, lats, lons);

            % create data array dimensions
            nrows = length(axes{2});
            ncols = length(axes{4});
            ntimes = length(refs);
            nlayers = length(output_layers);
            nfields = length(output_fields);

            % single-precision floats conserve memory and doesn't hurt accuracy
            % here because weather data is not accurate at the roundoff level anyway
            data = xarray(NaN(nrows, ncols, ntimes, nlayers, nfields, "single"), ...
                axes{:}, time = times, layer = output_layers, field = output_fields);

            for i_time = 1:ntimes
                inventory = refs(i_time).inventory;
                path = refs(i_time).local; 

                start = tic;
                [~, name, ext] = fileparts(path);
                fprintf("Reading from %%tempdir%%/%s (%s)... ", name + ext, refs(i_time).cycle); 
                for i_field = 1:nfields
                    field = output_fields(i_field);

                    filters.field = field;
                    filters.layer = output_layers;
                    filters.band = @isfinite;
                    messages = ncep.search(inventory, filters);

                    layers = inventory.layer(messages);

                    if length(layers) < length(output_layers)
                        warning("ncep:missing", "Missing %s data for %s", ...
                            field, mat2str(setdiff(output_layers, layers)));
                    end
                    
                    % Don't use <ismember>! If only a subset of the layers
                    % exist, the data can't be assigned to using ':'. This is
                    % the second time I've made this mistake.
                    [present, i_output, i_file] = intersect(output_layers, layers, "stable");
                    if isempty(present)
                        continue;
                    end

                    % read from file
                    field_data = readgeoraster(path, ...
                        Bands = inventory.band(messages), OutputType = "single");

                    % assign to output data
                    data{:, :, i_time, i_output, i_field} = ...
                        permute(field_data(data_indices{:}, i_file), [1 2 4 3 5]);
                end
                fprintf("done in %.1f sec\n", toc(start));
            end

            crs = info.CoordinateReferenceSystem;
        end

        function data = read_point(refs, lat, lon, params)
            arguments
                refs (1, :) ncep;
                lat (1,1) double;
                lon (1,1) double;
                params.field (1, :) pattern = wildcardPattern;
                params.layer (1, :) pattern = wildcardPattern;
            end

            [grid_data, crs] = refs.read(lats = lat + [-1 1], lons = lon + [-1 1], ...
                field = params.field, layer = params.layer);
            

            sample_values = ncep.latlon2coords(crs, lat, lon);
            sample_args = namedargs2cell(sample_values);
            data = grid_data.interp(sample_args{:});

            if all(ismember(["UGRD", "VGRD"], data.field)) && isa(crs, "projcrs")
                uv = data.align("field").pick{"field", ["UGRD", "VGRD"]};

                jacob = ncep.projjacob(crs, "geographic", lat, lon);
                jacob = jacob ./ vecnorm(jacob, 2, 1);

                winds = pagemtimes(jacob, uv);

                remaining_dims = repmat({':'}, 1, ndims(winds)-1);

                data.pick{"field", "UGRD"} = winds(1, remaining_dims{:});
                data.pick{"field", "VGRD"} = winds(2, remaining_dims{:});
            end
        end

        function delete(ref)
            if isfile(ref.local)
                delete(ref.local);
            end
        end
    end

    methods (Static, Access = public)
        % List valid models, products, cycles, depths given previous arguments
        % All arguments are optional, and the valid values for the next argument area always returned
        function [values, sources] = list(model, product, cycle, depth)
            % [values, sources] = list([model, product, cycle, depth])
            arguments
                model (1,1) string = missing;
                product (1,1) string = missing;
                cycle (1,1) duration = missing;
                depth (1,1) duration = missing;
            end
            
            err_id = "ncep:invalid_product";

            valid_models = ["gfs", "nam", "hrrr"];
            if ismissing(model)
                values = valid_models;
                sources = [];
                return;
            end

            switch model
                case "gfs"
                    valid_products = ["pgrb2.1p00", "pgrb2.0p50", "pgrb2.0p25"];
                    sources = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/" + ...
                        "gfs.%year$04d%month$02d%day$02d/%cycle$02d/atmos/" + ...
                        "gfs.t%cycle$02dz.%product$s.f%fcst$03d";
                case "nam"
                    valid_products = ["awphys", "conusnest.hiresf"];
                    sources = "https://noaa-nam-pds.s3.amazonaws.com/" + ...
                        "nam.%year$04d%month$02d%day$02d/" + ...
                        "nam.t%cycle$02dz.%product$s%fcst$02d.tm00.grib2";
                case "hrrr"
                    valid_products = "wrfprsf";
                    sources = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/" + ...
                        "hrrr.%year$04d%month$02d%day$02d/conus/" + ...
                        "hrrr.t%cycle$02dz.%product$s%fcst$02d.grib2";
                otherwise
                    error(err_id, "Unrecognized model '%s'", model);
            end

            % if no products specified, return valid options
            if ismissing(product)
                values = valid_products;
                return;
            end

            % validate product, then get valid cycles
            if ~ismember(product, valid_products)
                error(err_id, "Unsupported product '%s', must be one of \n%s", ...
                    product, ncep.printlist(valid_products));
            end

            switch model
                case "hrrr"
                    valid_cycles = hours(0:1:23);
                    hint = "hourly starting at 00Z";
                otherwise 
                    valid_cycles = hours(0:6:18);
                    hint = "6-hourly starting at 00Z";
            end

            % if no cycles specified, return valid options
            if ismissing(cycle)
                values = valid_cycles;
                return;
            end

            % validate cycle, then get valid depths
            if ~ismember(cycle, valid_cycles)
                error(err_id, "Invalid cycle %s: %s.%s is produced %s", ...
                    string(cycle), model, product, hint);
            end

            switch product
                case {"pgrb2.1p00", "pgrb2.0p50"}
                    valid_depths = hours(0:3:384);
                    hint = "3-hourly up to 384 hours";
                case "pgrb2.0p25"
                    valid_depths = [hours(0:1:120), hours(123:3:384)];
                    hint = "hourly up to 120 hours, 3-hourly up to 384 hours";
                case "awphys" 
                    valid_depths = [hours(0:1:36), hours(39:3:84)];
                    hint = "hourly up to 36 hours, 3-hourly up to 84 hours";
                case "conusnest.hiresf"
                    valid_depths = hours(0:1:60);
                    hint = "hourly up to 60 hours";
                case "wrfprsf"
                    if ismember(cycle, hours(0:6:18))
                        valid_depths = hours(0:1:48);
                        hint = "hourly up to 48 hours for 6-hourly cycles";
                    else
                        valid_depths = hours(0:1:18);
                        hint = "hourly up to 18 hours for hourly cycles";
                    end
                otherwise
                    error("This line should be unreachable");
            end

            % if no depth specified, return valid depths
            if ismissing(depth)
                values = valid_depths;
                return;
            end

            if ~ismember(depth, valid_depths)
                error(err_id, "Invalid depth %s: %s.%s is produced %s", ...
                    string(depth), model, product, hint);
            end
            
            values = depth;
        end

        function refs = forecast(model, product, validtime, reftime)
            arguments
                model (1,1) string;
                product (1,1) string;
                validtime (1,:) datetime {mustBeFinite, ncep.shouldHaveTimezone};
                reftime (1,1) datetime {mustBeFinite, ncep.shouldHaveTimezone};
            end

            reftime.TimeZone = "UTC";
            validtime.TimeZone = "UTC";
            if any(reftime > validtime)
                error("Reference time must occur on or before valid time");
            end
            ref_day = dateshift(reftime, "start", "day");

            candidate_cycles = ref_day + ncep.list(model, product);
            i_cycle = find(candidate_cycles <= reftime, 1, "last");
            nearest_cycle = candidate_cycles(i_cycle);

            assert(nearest_cycle <= reftime, "Internal error: " + ...
                "rounded-off cycle should occur on or before specified reference time");

            candidate_depths = ncep.list(model, product, nearest_cycle - ref_day);
            if isscalar(validtime)
                [~, i_depth] = min(abs(nearest_cycle + candidate_depths - validtime));
                depths = candidate_depths(i_depth);
            else
                i_fst = find(nearest_cycle + candidate_depths <= min(validtime), 1, "last");
                i_lst = find(nearest_cycle + candidate_depths >= max(validtime), 1, "first");
                depths = candidate_depths(i_fst:i_lst);
            end

            refs = arrayfun(@(dep) ncep(model, product, nearest_cycle, dep), depths);
        end

        function refs = analysis(model, product, validtime)
            arguments
                model (1,1) string;
                product (1,1) string;
                validtime (1,:) datetime {mustBeFinite, ncep.shouldHaveTimezone};
            end

            validtime.TimeZone = "UTC";

            valid_day = dateshift(validtime, "start", "day");
            candidate_days = min(valid_day):days(1):max(valid_day);
            cycles = ncep.list(model, product);
            candidate_datetimes = candidate_days + cycles(:); 
            % days is row, cyces is column - this becomes a matrix with all sums
            candidate_datetimes = sort(candidate_datetimes(:));

            if isscalar(validtime)
                [~, i_time] = min(abs(candidate_datetimes - validtime));
                datetimes = candidate_datetimes(i_time);
            else
                i_fst = find(candidate_datetimes <= min(validtime), 1, "last");
                i_lst = find(candidate_datetimes >= max(validtime), 1, "first");
                datetimes = candidate_datetimes(i_fst:i_lst);
            end
            
            refs = arrayfun(@(cyc) ncep(model, product, cyc, hours(0)), datetimes);
        end

        function jacob = projjacob(crs, mode, first, second)
            arguments
                crs projcrs;
                mode string;
                first double;
                second double;
            end

            if mode == "geographic"
                lat = first;
                lon = second;
                [x, y] = projfwd(crs, lat, lon);
            elseif mode == "planar"
                x = first;
                y = second;
                [lat, lon] = projinv(crs, x, y);
            end

            geo = crs.GeographicCRS;
            radius = mean([geo.Spheroid.SemimajorAxis, geo.Spheroid.SemiminorAxis]);
            delta = 1;
            jacob = zeros(2,2);
            [jacob(2,1), jacob(1,1)] = projinv(crs, x + delta, y);
            [jacob(2,2), jacob(1,2)] = projinv(crs, x, y + delta);
            jacob = jacob - [lon; lat];
            jacob = deg2rad(jacob) .* radius .* [cosd(lat); 1];
        end

        function [rows] = search(index, filters)
            arguments (Input)
                index (:, :) table;
                filters (1,1) struct;
            end
            arguments (Output)
                rows (:,1) double {mustBeInteger, mustBePositive}
            end

            err_id = "weathergrid:invalidFilter";
            fields = string(fieldnames(filters));
            notpresent = setdiff(fields, index.Properties.VariableNames);
            if ~isempty(notpresent)
                % NOTE: may be better to emit error here
                warning(err_id, "Ignoring fields %s: not present in index", mat2str(notpresent));
            end
            fields = setdiff(fields, notpresent);

            rows = true(height(index), 1);

            for i_field = 1:length(fields)
                field = fields(i_field);
                filter = filters.(field);
                if isa(filter, "function_handle")
                    result = filter(index.(field));
                    if ~islogical(result)
                        error(err_id, "Filter <%s> returned outputs of type %s (%s expected)", ...
                            func2str(filter), class(result), "logical");
                    end
                    if length(result) ~= length(rows)
                        error(err_id, "Filter <%s> on field %s returned %d outputs (%d expected)", ...
                            func2str(filter), field, numel(result), numel(rows))
                    end
                else
                    result = matches(index.(field), filter);
                end

                rows = rows & result;
            end
            rows = find(rows);
        end
    end

    methods (Access = protected)
        function download(ref, opts)
            arguments
                ref (1,1) ncep
                opts (1,1) struct;
            end


            if ismissing(ref.local)
                % in theory, a specific folder is unnecessary becasue delete() gets
                % rid of the local file, but just in case, all the .grib files
                % should get shoved into this temp folder so they can be deleted
                % all at once if that fails
                tn = tempname;
                [folder, name, ~] = fileparts(tn);
                download_folder = fullfile(folder, "grib", name);
                [status, msg] = mkdir(download_folder);
                if ~status
                    error("%d: Unable to create temporary directory for downloading data (%s)", ...
                        status, msg);
                end

                [~, name, ext] = fileparts(ref.origin);
                filename = name + ext;
                ref.local = fullfile(download_folder, filename);
            end

            % find GRIB messages to download
            opts.band = @isnan;
            messages = ncep.search(ref.inventory, opts);
            if isempty(messages)
                fprintf("No messages to download\n");
                return;
            end

            % status output

            [~, name, ext] = fileparts(ref.origin);
            fprintf("Located %d messages on NOMADS at %s\n", length(messages), name + ext);
            data_file_id = fopen(ref.local, "a");
            if data_file_id < 0
                error("ncep:cantwrite", "Unable to open local file for writing");
            end
            oc = onCleanup(@() fclose(data_file_id));

            % download messages
            t_download = tic;
            char_previous = 0;
            bytes_downloaded = 0;
            for i_msg = 1:length(messages)
                t_message = tic;

                % range header start-end or start- depending on whether the message is the last one
                message = messages(i_msg);
                if message == height(ref.inventory)
                    opts = ncep.range_header(ref.inventory.start(message));
                elseif ref.inventory.start(message) == ref.inventory.start(message+1);
                    % this is a submessage, skip to the last one
                    continue;
                else
                    opts = ncep.range_header(ref.inventory.start(message), ...
                        ref.inventory.start(message+1)-1);
                end

                try
                    % somehow about twice as fast as remote-file fopen() and fread()???
                    bytes = webread(ref.origin, opts);
                    fwrite(data_file_id, bytes);
                    bytes_downloaded = bytes_downloaded + length(bytes);
                catch mex
                    if mex.identifier == "MATLAB:webservices:HTTP404StatusCodeError"
                        error("ncep:notfound", "Not found (404): %s", url);
                    else
                        rethrow(mex);
                    end
                end

                % status output
                [num, unit] = ncep.format_filesize(length(bytes));
                status_message = sprintf("Downloaded message %3d of %3d (%.1f %s) in %.3g sec\n", ...
                    i_msg, length(messages), num, unit, toc(t_message));

                fprintf(repmat('\b', 1, char_previous) + status_message);
                char_previous = strlength(status_message);
            end

            % register addition of messages as bands
            % using (band) becasue it is what readgeoraster() uses
            ref.inventory.band(messages) = 1:length(messages) + max([0; ref.inventory.band]);
            [num, unit] = ncep.format_filesize(bytes_downloaded);
            fprintf("Downloaded %d messages (%.1f %s) in %.3g sec\n", ...
                length(messages), num, unit, toc(t_download));
        end
    end

    methods (Static, Access = protected) 
        % Read NCEP index from url
        function tbl = read_ncep_index(url)
            arguments
                url (1,1) string;
            end

            file_columns = ["message", "start", "cycle", "field", "layer", "fcst", "attr"];
            opts = delimitedTextImportOptions(NumVariables = length(file_columns), ...
                Delimiter = ":", VariableNames = file_columns);
            opts = setvaropts(opts, "message", Type = "double");
            opts = setvaropts(opts, "start", Type = "int32");
            opts = setvaropts(opts, "cycle", Type = "datetime", ...
                InputFormat = "'d='yyyyMMddHH", TimeZone = "UTC", ...
                DatetimeFormat = "yyyy-MM-dd HHz");
            opts = setvaropts(opts, ["field", "layer", "fcst", "attr"], Type = "string");

            tbl = readtable(url, opts);
        end

        function tbl = patch_index(tbl)
            band = NaN(height(tbl), 1);
            tbl = addvars(tbl, band, After = "message");
        end

        % Search table for filters specified in structure

        % Geographically crop planar or geographic raster reference, with
        % sane behavior for out-of-range or Infinite limits. Return the indices 
        % into the original raster so that this does not need to be repeatedly re-calculated.
        function [raster, indices, axes] = crop_raster(raster, lats, lons)
            % [raster, indices, axes] = crop_raster(raster, lats, lons)
            % Inputs
            %   raster  (Cell or Posting reference)  raster to crop
            %   lats    (double)    [min max] latitude limits, in any order (+/- Inf are valid)
            %   lons    (double)    [min max] longitude limits, in any order (+/- Inf are valid)
            %
            % Outputs
            %   raster  (Cell or Posting reference)  cropped raster reference
            %   indices {rows-logical, cols-logical} 
            %               Cell of logical vectors that will index into the data
            %   axes    {rows-name, rows-values, cols-name, cols-values} 
            %               Cell of axes to pass into <xarray> defining names
            %               and values of the cropped coordinate axes
            arguments (Input)
                raster                 
                lats (1,2) double;
                lons (1,2) double;
            end
            arguments (Output)
                raster;
                indices (1,2) cell;
                axes (1,4) cell;
            end
            raster_type = raster.CoordinateSystemType;
            
            % retreiving the grid (vectors) for each raster and comparing them
            % to lats/lons to get the inputs to (map/geo)crop lets this
            % function more succinctly support infinite bounds for lats or lons
            % if required
            switch raster_type
                case "planar"
                    [x, y] = worldGrid(raster, "gridvectors");

                    % The lat/lon extrema can happen anywhere in the raster,
                    % depending on the projection. So, we project the entire
                    % coordinate grid to geographic coordinates and find the
                    % extrema in that. 

                    [x_grid, y_grid] = worldGrid(raster); 
                    [lat, lon] = projinv(raster.ProjectedCRS, x_grid, y_grid);

                    % filter x/y grids based on lat/lon limits
                    x_filter = (min(lons) <= lon & lon <= max(lons));
                    y_filter = (min(lats) <= lat & lat <= max(lats));
                    x_in_range = x_grid(x_filter & y_filter);
                    y_in_range = y_grid(x_filter & y_filter);

                    % calculate x/y limits
                    x_lims = [min(x_in_range, [], "all") max(x_in_range, [], "all")];
                    y_lims = [min(y_in_range, [], "all") max(y_in_range, [], "all")];
                    if isempty(x_lims) || isempty(y_lims)
                        error("No points found in specified range");
                    end

                    % mapcrop
                    [~, raster] = mapcrop(zeros(length(y), length(x)), raster, x_lims, y_lims);

                    % recalculate limits because mapcrop may leave an extra element on either side
                    row_idx = raster.YWorldLimits(1) <= y & y <= raster.YWorldLimits(2);
                    col_idx = raster.XWorldLimits(1) <= x & x <= raster.XWorldLimits(2);
                        
                    % return indices 
                    indices = {row_idx, col_idx};
                    axes = {"y", y(row_idx), "x", x(col_idx)};
                case "geographic"
                    % lat/lon filtering of geographic coordinates is much more straightforward
                    [lat, lon] = geographicGrid(raster, "gridvectors");

                    lat_in_range = lat(min(lats) <= lat & lat <= max(lats));
                    lon_in_range = lon(min(lons) <= lon & lon <= max(lons));
                    lon_lims = [min(lon_in_range, [], "all") max(lon_in_range, [], "all")];
                    lat_lims = [min(lat_in_range, [], "all") max(lat_in_range, [], "all")];

                    [~, raster] = geocrop(zeros(length(lat), length(lon)), raster, lat_lims, lon_lims);

                    % same "extra element" patch as before
                    row_idx = raster.LatitudeLimits(1) <= lat & lat <= raster.LatitudeLimits(2);
                    col_idx = raster.LongitudeLimits(1) <= lon & lon <= raster.LongitudeLimits(2);

                    indices = {row_idx, col_idx};
                    axes = {"lat", lat(row_idx), "lon", lon(col_idx)};
                otherwise
                    error("Unrecognized coordinate system type '%s'", raster_type);
            end
        end

        function values = latlon2coords(crs, lat, lon)
            if isa(crs, "geocrs")
                values.lat = lat;
                values.lon = lon;
            elseif isa(crs, "projcrs")
                [x, y] = projfwd(crs, lat, lon);
                values.y = y;
                values.x = x;
            else
                error("Invalid coordinate reference system");
            end
        end


        % Calculate engineering-unit file sizes (k/M/G/T...)
        function [nums, units] = format_filesize(sizes)
            % [nums, units] = format_filesize(sizes)
            arguments (Input)
                sizes (1,:) double;
            end
            arguments (Output)
                nums (1,:) double;
                units (1,:) string;
            end

            pwr = floor(log10(sizes)/3)*3;
            units_map = dictionary([-Inf 0 3 6 9 12], ["B", "B", "kB", "MB", "GB", "TB"]);
            units = units_map(pwr);
            nums = sizes./10.^pwr;
            nums(~isfinite(nums)) = 0;
        end

        % Adds name-based indexing to sprintf's %<number>$<spec>
        % capability, so that %<name>$<spec> will use the corresponding
        % Name-Value pair for formatting.
        function [str] = sprintnf(str, name, rep)
            % str = sprintnf(str, Name, Value, ...)
            % 
            % Example:
            % sprintnf("%year$4d-%month$s-%day$2d", ...
            %   month = "Jan", year = 1970, day = 1);

            arguments
                str (:, 1) string;
            end
            arguments (Repeating)
                name (1, 1) string;
                rep (:, 1);
            end
            
            % match c-style format specifier
            formatop = optionalPattern("-" | "+" | " " | "0" | "#") + ...
                optionalPattern(digitsPattern) + optionalPattern(".") + ...
                optionalPattern(digitsPattern) + lettersPattern(1);

            % match single % followed by name, followed by spec
            repfield = lookAheadBoundary("%") + "%" + lettersPattern + ...
                "$" + wildcardPattern + formatop;

            % convert strings to numbers based on order of [name]
            locs = extract(str, repfield);
            strnames = extractBetween(locs, "%", "$");
            [strnames, ~, indices] = intersect(strnames, [name{:}], "stable");
            str = replace(str, strnames, string(indices));

            % ordinary sprintf
            for i_str = 1:numel(str)
                str(i_str) = sprintf(str, rep{:}); 
            end
        end

        % Get the HTTP header for a range-request
        function wopts = range_header(start, fin)
            % wopts = range_header(start[, fin])
            arguments
                start (1,1) double;
                fin (1,1) double = Inf;
            end
            if isfinite(fin)
                range_text = sprintf("bytes=%d-%d", start, fin);
            else
                range_text = sprintf("bytes=%d-", start);
            end
            wopts = weboptions(RequestMethod = "get", ...
                HeaderFields = ["Range", range_text]);
        end

        % Validate that table has specified variables
        function mustHaveCols(tbl, vars)
            nocol = setdiff(vars, tbl.Properties.VariableNames);
            tblname = inputname(1);
            if tblname == ""
                tblname = "table";
            end
            if ~isempty(nocol)
                mex = MException("ncep:corrupt", "%s is missing columns %s", ...
                    tblname, mat2str(nocol));
                throwAsCaller(mex);
            end
        end

        function shouldHaveTimezone(time)
            if isempty(time.TimeZone)
                warning("ncep:timezone", ...
                    "datetime input has no assigned time zone, defaulting to UTC");
            end
        end


        % convert (names, values) Repeating arguments into a structure
        function stru = namevalue2struct(names, values)
            for i_field = 1:length(names)
                stru.(names{i_field}) = values{i_field};
            end
        end

        function [stru, val] = strip(stru, field, deflt)
            if isfield(stru, field)
                val = stru.(field);
                stru = rmfield(stru, field);
            elseif nargin == 2
                error("Struct has no field '%s' and no default value specified", field);
            elseif nargin == 3
                val = deflt;
            end
        end

        function str = printlist(strs, fmt)
            arguments
                strs (:,1) string;
                fmt (1,1) string = "\t""%s""\n"
            end

            str = compose(fmt, strs).join("");
        end

    end
end
