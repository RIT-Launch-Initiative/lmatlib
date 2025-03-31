%% Retrieve and read NCEP GRIB2 files
% Uses the AWS S3 bucket maintained by the Open Data Dissemination Program
% BRIEF
%   ncep.list           list models and products
%   ncep(...)           construct single reference from exact cycle and forecast
%   ncep.fcst(...)      construct array of forecasts from reference time and range of valid-times
%   ncep.anl(...)       construct array of analyses from range of valid-times
% 
% 

classdef ncep < handle
    properties (GetAccess = public, SetAccess = protected)
        model (1,1) string = string(missing);
        product (1,1) string = string(missing);
        cycle (1,1) datetime = datetime(missing);
        forecast (1,1) duration = duration(missing);
        member (1,1) string = string(missing);
    end

    properties (Access = public)
        filename (1,1) string = string(missing);
        local_path (1,1) string = string(missing);
        remote_path (1,1) string = string(missing);
        cache; 
    end
    
    properties (Constant, Access = protected)
        defs dictionary = ncep.data_definitions;
    end

%% BASIC I/O
    methods
        function obj = ncep(model, product, cycle, forecast)
            % ref = ncep(model, product, cycle, forecast[, roundtime])
            %   model   (string)    NWP model 
            %   grid    (string)    Data product
            %   cycle   (datetime)  Model output cycle (usually hourly or 6-hourly, 00 UTC)
            %   forecast(duration)  Forecast extent (usually hourly or 3-hourly)
            arguments
                model (1,1) string;
                product (1,1) string;
                cycle (1,1) datetime {mustBeFinite};
                forecast (1,1) duration {mustBeFinite, mustBeNonempty};
            end
            
            [modeldef, productdef] = ncep.validate_product(model, product);
            obj.model = model;
            obj.product = product;

            % Validate cycle --------------------------------------------------
            if isempty(cycle.TimeZone)
                warning("ncep:timezone", ...
                    "Model cycle does not have assigned time zone. Defaulting to UTC.");
            end

            cycle.TimeZone = "UTC";
            cycle_day = dateshift(cycle, "start", "day");
            cycle_time = cycle - cycle_day;
            valid_cycles = productdef{"cycles"};

            if ~ismember(cycle_time, valid_cycles)
                error("ncep:invalid", "Invalid model cycle %s: %s %s is produced %s", ...
                    cycle, model, product, productdef{"cycles_hint"});
            end

            cycle = cycle_day + cycle_time;
            cycle.Format = "dd-MMM-yyyy HH z";
            obj.cycle = cycle;

            % Validate forecasts ----------------------------------------------
            valid_fcst = productdef{"forecasts"};
            invalid = setdiff(forecast, valid_fcst);
            if ~isempty(invalid)
                error("Invalid model forecasts(s) %s: %s %s is produced %s", ...
                    mat2str(string(invalid)), obj.model, obj.product, ...
                    productdef{"forecasts_hint"});
            end
            obj.forecast = forecast;

            kwargs = {"product", obj.product, "cycle", obj.cycle.Hour, "fcst", hours(obj.forecast), ...
                "year", obj.cycle.Year, "month", obj.cycle.Month, "day", obj.cycle.Day};

            obj.filename = ncep.sprintnf(modeldef{"filename"}, kwargs{:});
            obj.remote_path = fullfile(ncep.sprintnf(modeldef{"url"}, kwargs{:}), obj.filename);

            obj.local_path = tempname;
            obj.cache = matfile(obj.local_path + ".mat", Writable = true);
        end

        function attach(refs, folder)
            arguments
                refs (1, :) ncep {mustBeNonempty};
                folder (1,1) string {mustBeFolder};
            end

            for i_ref = 1:length(refs)
                refs(i_ref).set_local_folder(folder);
            end
        end

        function download(refs, names, values)
            % Download subset of messages to local storage
            % refs.download(folder, Name, Value)
            %   folder  (string)    parent folder to download to
            %                       subfolders yyyy-MM-dd are created
            %   Name    (string)    inventory column to filter by
            %                       usually "field" and "layer"
            %   Value   (...)       filtering value
            %                       String or pattern array for text columns
            %                       Exact matches for all others (numeric, datetime)
            arguments
                refs (1, :) ncep;
            end
            arguments (Repeating)
                names (1,1) string;
                values (:, 1);
            end

            for i_ref = 1:length(refs)
                ref = refs(i_ref);
                inventory = ref.inventory;
                
                stacked = [names; values];
                indices = ncep.search(inventory, stacked{:});
                messages = inventory.message(indices);

                msg_to_download = ismember(inventory.message, messages) & isnan(inventory.band);

                n_to_download = sum(msg_to_download);
                fprintf("Located %d messages at %s\n", length(messages) - n_to_download, ref.local_path);
                if n_to_download == 0
                    continue;
                end
                info = ncep.add_messages(ref.local_path, ref.remote_path, inventory(msg_to_download, :));

                % register addtion of new bands to file
                added_bands = max([0; inventory.band]) + (1:n_to_download);
                inventory.band(msg_to_download) = added_bands;

                if max(inventory.band) ~= height(info.Metadata)
                    delete(ref.local_path);
                    delete(ref.cache.Proeprties.Source);
                    error("Inventory and RasterInfo disagree on available file bands for %s", ...
                        ref.local_path);
                end

                ref.cache.inventory = inventory;
            end
        end
        
        function [data, crs] = read(obj, params)
            arguments
                obj (1, :) ncep;
                params.lats (1,2) double = [-Inf Inf];
                params.lons (1,2) double = [-Inf Inf];
                params.field (1, :) pattern = wildcardPattern;
                params.layer (1, :) pattern = wildcardPattern;
            end

            template = obj.inventory;
            template = template(isfinite(template.band), :);
            msg_indices = ncep.search(template, ...
                field = params.field, layer = params.layer);
            output_fields = unique(template.field(msg_indices));
            output_layers = unique(template.layer(msg_indices));
            times = [obj.cycle] + [obj.forecast];

            info = georasterinfo(obj(1).local_path);
            [~, data_indices, axes] = ncep.crop_raster(...
                info.RasterReference, ...
                params.lats, params.lons);

            nrows = length(axes{2});
            ncols = length(axes{4});
            ntimes = length(obj);
            nlayers = length(output_layers);
            nfields = length(output_fields);

            data = xarray(NaN(nrows, ncols, ntimes, nlayers, nfields, "single"), ...
                axes{:}, time = times, layer = output_layers, field = output_fields);

            for i_time = 1:ntimes
                inventory = obj(i_time).inventory;
                inventory = inventory(isfinite(inventory.band), :);
                path = obj(i_time).local_path; 

                start = tic;
                fprintf("Reading from %s... ", path); 
                for i_field = 1:nfields

                    field = output_fields(i_field);
                    msg_indices = ncep.search(inventory, ...
                        field = field, layer = output_layers);
                    layers = inventory.layer(msg_indices);

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

                    % read and assign to output
                    field_data = readgeoraster(path, ...
                        Bands = inventory.band(msg_indices), OutputType = "single");
                    data{:, :, i_time, i_output, i_field} = ...
                        permute(field_data(data_indices{:}, i_file), [1 2 4 3 5]); 
                end
                fprintf("done in %.1f sec\n", toc(start));
            end


            crs = info.CoordinateReferenceSystem;
        end
        
        function varargout = inventory(nc)
            % Retrieve inventor(ies) from local or remote resource
            arguments
                nc (1,:) ncep {mustBeNonempty};
            end
            
            nout = max(1, nargout);
            varargout = cell(1, nout);
            for i_nc = 1:nout
                if isfile(nc(i_nc).cache.Properties.Source)
                    varargout{i_nc} = nc(i_nc).cache.inventory;
                else
                    varargout{i_nc} = ncep.read_index(nc(i_nc).remote_path);  
                end
            end
        end

        function varargout = raster(nc)
            nout = max(1, nargout);
            varargout = cell(1, nout);
            for i_nc = 1:nout
                if isfile(nc(i_nc).local_path)
                    varargout{i_nc} = nc(i_nc).cache.raster;
                else
                    error("File not downloaded, no raster information available") 
                end
            end
            
        end
    end

    methods (Static)
        function list
            models = ncep.defs.keys;
            for i_m = 1:length(models)
                mdl = models(i_m);
                fprintf("%s: %s\n", mdl, ncep.defs{mdl}{"name"});
                pdefs = ncep.defs{mdl}{"products"};
                products = pdefs.keys;
                for i_p = 1:length(products)
                    pd = products(i_p);
                    fprintf("\t%s: %s\n", pd, pdefs{pd}{"description"});
                end
                fprintf("\n");
            end
        end

        function refs = fcst(model, product, reftime, validrange)
            % Create a forecast for any reference time and range
            % refs = ncep.series(model, product, reftime, validrange)
            %   model   (string)        same as ncep() constructor (see ncep.list)
            %   product (string)        same as ncep() constructor (see ncep.list)
            %   reftime (datetime)      approximate reference time
            %                           uses the closest model cycle on or before this time
            %                           (e.g. 2024-01-01 05:30Z pushed back to 2024-01-01 00:00Z)
            %                           Accounts for time zone difference, assumes UTC if not set
            %   validrange (datetime)   (1x2) required forecast time range (min, max)
            %                           uses the narrowest range of model forecasts surrounding this range
            %                           (e.g 2024-01-01 [05:30Z 14:00Z] gives
            %                           2024-01-01 [03 06 09 12 15]Z) for a 3-hourly model
            %                           Accounts for time zone difference, assumes UTC if not set
            %   
            %   OUTPUT
            %   refs  Array of ncep references
            arguments
                model (1,1) string; 
                product (1,1) string; 
                reftime (1,1) datetime {mustBeFinite}; 
                validrange (1,2) datetime {mustBeFinite};
            end

            [~, productdef] = ncep.validate_product(model, product);

            reftime.TimeZone = "UTC";
            validrange.TimeZone = "UTC";

            cycle_day = dateshift(reftime, "start", "day");
            cycle_time = reftime - cycle_day;
            valid_cycles = productdef{"cycles"};
            i_rounded = find(cycle_time >= valid_cycles, 1, "last");
            cycle_time = valid_cycles(i_rounded);
            cycle = cycle_day + cycle_time;

            assert(cycle <= reftime, ...
                "Model cycle must occur at or before specified reference time");

            fcst_range = validrange - cycle;
            valid_fcst = productdef{"forecasts"};
            i_start = find(fcst_range(1) >= valid_fcst, 1, "last");
            i_end = find(fcst_range(2) <= valid_fcst, 1, "first");
            forecasts = valid_fcst(i_start:i_end);
            
            assert((cycle + forecasts(1)) <= validrange(1), ...
                "Forecast series must start at or before first valid time");
            assert((cycle + forecasts(end)) >= validrange(2), ...
                "Forecast series must end at or after second valid time");

            refs = arrayfun(@(fcst) ncep(model, product, cycle, fcst), forecasts);
        end

        function refs = anl(model, product, valid_range)
            arguments
                model (1,1) string; 
                product (1,1) string; 
                valid_range (1,2) datetime {mustBeFinite};
            end

            [~, productdef] = ncep.validate_product(model, product);
            product_cycles = productdef{"cycles"};

            % round times up/down to create boundaries
            valid_range.TimeZone = "UTC";
            boundary_days = dateshift(valid_range, "start", "day");

            cand_start_times = boundary_days(1) + product_cycles; % candidate start-times
            range_start = cand_start_times(find(valid_range(1) >= cand_start_times, 1, "last"));
            cand_end_times = boundary_days(2) + product_cycles(:) + [days(0), days(1)]; % candidate end-times
            % must increment cycle along column, day along row, becasue find() index is column-major
            range_end = cand_end_times(find(valid_range(2) <= cand_end_times, 1, "first"));

            assert((range_start <= valid_range(1)) && (valid_range(2) <= range_end), ...
                "Analysis series must envelope specified valid range");
            
            % candidate times
            range_days = dateshift([range_start range_end], "start", "day");
            cand_times = (range_days(1):days(1):range_days(2)) + product_cycles';
            selected_cycles = cand_times((range_start <= cand_times) & (cand_times <= range_end));

            refs = arrayfun(@(cycle) ncep(model, product, cycle, hours(0)), selected_cycles);
        end
    end

    methods (Access = protected)
        function set_local_folder(ref, folder)
            arguments
                ref (1,1) ncep;
                folder (1,1) {mustBeFolder};
            end

            dest = fullfile(folder, string(ref.cycle, "yyyy-MM-dd"));
            [status, ~, ~] = mkdir(dest);
            if ~status
                error("Unable to create destination folder %s", dest);
            end
            
            ref.local_path = fullfile(dest, ref.filename);
            ref.cache = matfile(ref.local_path + ".mat", Writable = true);
        end
    end

%% DATA IMPORT PRIMITIVES
    methods (Static, Access = public)

        function [info] = add_messages(file, url, index)
            % info = add_messages(url, index)
            % Download messages specified in <index> to temporary file
            % Return RasterInfo

            arguments
                file (1,1) string;
                url (1,1) string;
                index table {ncep.mustHaveCols(index, ["message", "start", "len"]), ...
                    mustBeNonempty};
            end

            t_start = tic;

            % TODO validate that all components of a submessage have been specified

            flocal = fopen(file, "a");
            if flocal < 0
                error("Unable to open temporary file for writing: %d", flocal);
            end
            oc = onCleanup(@() fclose(flocal)); % close file no matter what

            
            fprintf("Downloading %d messages (%s)...", ...
                height(index), ncep.format_filesize(sum(index.len)));

            for i_msg = 1:height(index)
                start = index.start(i_msg);
                len = index.len(i_msg);
                if len > 0
                    range_text = sprintf("bytes=%d-%d", start, start+len-1);
                    opts = weboptions(RequestMethod = "get", ...
                        HeaderFields = ["Range", range_text], ...
                        Timeout = 30);
                    try
                        % somehow about twice as fast as remote-file fopen() and fread()???
                        bytes = webread(url, opts);
                    catch mex
                        if mex.identifier == "MATLAB:webservices:HTTP404StatusCodeError"
                            error("ncep:notfound", "Not found (404): %s", url);
                        else
                            rethrow(mex);
                        end
                    end

                    if length(bytes) ~= len
                        error("ncep:corrupt", "Expected %d bytes, read %d", len, bytes);
                    end
                    fwrite(flocal, bytes);
                end
            end


            info = georasterinfo(file);

            time = toc(t_start);
            fprintf(" finished in %.1f sec\n", time);
        end

        function [tbl] = read_index(grib_url)
            % Read NCEP .grib file inventory (.idx) from remote URL
            % Populate default values for addtional arguments
            arguments
                grib_url (1,1) string;
            end

            start = tic; 
            index_url = grib_url + ".idx";
            fprintf("Reading index from %s...", grib_url + ".idx");
            
            % See https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/default_inv.html
            file_columns = ["message", "start", "cycle", "field", "layer", "fcst", "attr"];
            opts = delimitedTextImportOptions(NumVariables = length(file_columns), ...
                Delimiter = ":", VariableNames = file_columns);
            opts = setvaropts(opts, "message", Type = "double");
            opts = setvaropts(opts, "start", Type = "int32");
            opts = setvaropts(opts, "cycle", Type = "datetime", ...
                InputFormat = "'d='yyyyMMddHH", TimeZone = "UTC", ...
                DatetimeFormat = "yyyy-MM-dd HHz");
            opts = setvaropts(opts, ["field", "layer", "fcst", "attr"], Type = "string");

            tbl = readtable(index_url, opts);

            % add length column, for convenience
            % used by sync_bands for range-requesting
            info = dir(grib_url);
            len = [tbl.start(2:end); info.bytes] - tbl.start;
            assert(all(len >= 0), "All message lengths must be nonnegative");
            tbl = addvars(tbl, len, After = "start");

            band = NaN(height(tbl), 1);
            tbl = addvars(tbl, band, After = "message");

            tbl.Properties.Description = "Index for NCEP GRIB2 file at " + grib_url;

            time = toc(start);
            fprintf(" done in %.1f sec\n", time);
        end

        function [indices] = search(index, names, values)
            arguments
                index table
            end
            arguments (Repeating)
                names (1,1) string;
                values (:, 1);
            end

            mask = true(height(index), 1);

            for i_filter = 1:length(names)
                col = names{i_filter};

                if index.Properties.VariableTypes(col) == "string"
                    % for string-type columns, use pattern-matching
                    filter = matches(index.(col), values{i_filter});
                else 
                    % otherwise use <ismember> for exact matching
                    filter = ismember(index.(col), values{i_filter});
                end

                mask = mask & filter;
            end
            
            indices = find(mask);
        end
    end

%% MISCALLANEOUS
    methods (Static, Access = public)
        function [modeldef, productdef] = validate_product(model, product)
            if ~ismember(model, ncep.defs.keys)
                mex = MException("ncep:invalid", ...
                    "Model ""%s"" not recognized. Models supported:\n\t%s", ...
                    model, ncep.defs.keys.join(", "));
                throwAsCaller(mex);
            end
            modeldef = ncep.defs{model};
            if ~ismember(product, modeldef{"products"}.keys)
                mex = MException("ncep:invalid", ...
                    """%s"" not recognized. %s supports\n\t%s", ...
                    product, model, modeldef{"products"}.keys.join(", "));
                throwAsCaller(mex);
            end
            productdef = modeldef{"products"}{product};
        end

        function [defs] = data_definitions
            % Create data definitions
            % Nested dictionary defining models and products

            model = configureDictionary("string", "cell");
            model{"filename"} = string(missing);
            model{"url"} = string(missing);

            product = configureDictionary("string", "cell");
            product{"description"} = string(missing);
            product{"cycles"} = hours(0:6:18);
            product{"cycles_hint"} = "6-hourly starting at 00Z";
            product{"forecasts"} = NaT;
            product{"forecasts_hint"} = string(missing);
            product{"members"} = string([]);
            product{"members_hint"} = string(missing);

            onedeg = product;
            onedeg{"description"} = "GFS 1.00 deg lat/lon commonly-used atmospheric fields";
            onedeg{"forecasts"} = hours(0:3:384);
            onedeg{"forecasts_hint"} = "3-hourly up to 384 hours";

            halfdeg = onedeg;
            halfdeg{"description"} = "GFS 0.50 deg lat/lon commonly-used atmospheric fields";

            quartdeg = onedeg;
            quartdeg{"description"} = "GFS 0.25 deg lat/lon commonly-used atmospheric fields";
            quartdeg{"forecasts"} = [hours(0:1:120), hours(123:3:384)];
            quartdeg{"forecasts_hint"} = "hourly up to 120 hours, 3-hourly up to 384 hours";

            gfs = model;
            gfs{"products"} = dictionary(["pgrb2.1p00", "pgrb2.0p50", "pgrb2.0p25"], ...
                {onedeg halfdeg quartdeg});
            gfs{"name"} = "Global Forecast System";
            gfs{"filename"} = "gfs.t%cycle$02dz.%product$s.f%fcst$03d";
            % gfs{"url"} = "s3://noaa-gfs-bdp-pds/" + ...
            gfs{"url"} = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/" + ...
                "gfs.%year$04d%month$02d%day$02d/%cycle$02d/atmos";

            ehalfdeg = halfdeg;
            ehalfdeg{"description"} = "GEFS 0.50 deg lat/lon commonly-used atmospheric fields";
            ehalfdeg{"members"} = ["c00", "avg", "spr", compose("p%02d", (1:30)')'];
            ehalfdeg{"members_hint"} = "control: c00, mean: avg, spread: spr, member: pNN";

            gefs = model;
            gefs{"products"} = dictionary("pgrb2a.0p50", {ehalfdeg});
            gefs{"name"} = "Global Ensemble Forecast System";
            gefs{"filename"} = "ge%member$s.t%cycle$02dz.%product$s.f%fcst$03d";
            gefs{"url"} = "https://noaa-gefs-pds.s3.amazonaws.com/" + ...
                "gefs.%year$04d%month$02d%day$02d/%cycle$02d/atmos/pgrb2ap5";

            % NAM definitions
            twelve = product;
            twelve{"description"} = "NAM 12-km Lambert (CONUS) pressure fields";
            twelve{"forecasts"} = [hours(0:1:36), hours(39:3:84)];
            twelve{"forecasts_hint"} = "hourly up to 36 hours, 3-hourly up to 84 hours";

            three = product;
            three{"description"} = "NAM 3-km Lambert (CONUS high-resolution) pressure fields";
            three{"forecasts"} = hours(0:1:60);
            three{"forecasts_hint"} = "hourly up to 60 hours";

            nam = model;
            nam{"products"} = dictionary(...
                ["awphys", "conusnest.hiresf"], ...
                {twelve, three});
            nam{"name"} = "North American Mesoscale";
            nam{"filename"} = "nam.t%cycle$02dz.%product$s%fcst$02d.tm00.grib2";
            % nam{"url"} = "s3://noaa-nam-pds/" + ...
            nam{"url"} = "https://noaa-nam-pds.s3.amazonaws.com/" + ...
                "nam.%year$04d%month$02d%day$02d";

            % HRRR definitions
            hrrr_product = product;
            hrrr_product{"description"} = "HRRR 3-km Lambert (CONUS high-resolution) pressure fields";
            hrrr_product{"forecasts"} = hours(0:1:48);
            hrrr_product{"forecasts_hint"} = "hourly up to 48 hours";
            hrrr_product{"cycles"} = hours(0:1:23);
            hrrr_product{"cycles_hint"} = "hourly";

            hrrr = model;
            hrrr{"products"} = dictionary("wrfprsf", {hrrr_product}); % USER FACING NAMES
            hrrr{"name"} = "High-Resolution Rapid Refresh";
            hrrr{"filename"} = "hrrr.t%cycle$02dz.%product$s%fcst$02d.grib2";
            % hrrr{"url"} = "s3://noaa-hrrr-bdp-pds/" + ...
            hrrr{"url"} = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/" + ...
                "hrrr.%year$04d%month$02d%day$02d/conus";

            % USER FACING NAMES
            defs = dictionary(...
                ["gfs", "gefs", "nam", "hrrr"], ...
                {gfs, gefs, nam, hrrr}); 
        end

        function [raster, indices, axes] = crop_raster(raster, lats, lons)
            arguments
                raster 
                lats (1,2) double;
                lons (1,2) double;
            end
            raster_type = raster.CoordinateSystemType;
            
            % retreiving the grid (vectors) for each raster and comparing them
            % to lats/lons to get the inputs to (map/geo)crop lets this
            % function more succinctly support infinite bounds for lats or lons
            % if required
            switch raster_type
                case "planar"
                    [x, y] = worldGrid(raster, "gridvectors");

                    % the lat/lon extrema can happen anywhere in the raster, depending on the projection
                    % so we project and compare based on the entire grid
                    [x_grid, y_grid] = worldGrid(raster); 
                    [lat, lon] = projinv(raster.ProjectedCRS, x_grid, y_grid);
                    x_in_range = x_grid(min(lons) <= lon & lon <= max(lons));
                    y_in_range = y_grid(min(lats) <= lat & lat <= max(lats));
                    x_lims = [min(x_in_range, [], "all") max(x_in_range, [], "all")];
                    y_lims = [min(y_in_range, [], "all") max(y_in_range, [], "all")];

                    [~, raster] = mapcrop(zeros(length(y), length(x)), raster, x_lims, y_lims);

                    % need to use the new limits instead of x_lims/y_lims because mapcrop
                    % leaves an extra element on either side
                    row_idx = raster.YWorldLimits(1) <= y & y <= raster.YWorldLimits(2);
                    col_idx = raster.XWorldLimits(1) <= x & x <= raster.XWorldLimits(2);
                        
                    indices = {row_idx, col_idx};
                    axes = {"y", y(row_idx), "x", x(col_idx)};
                case "geographic"
                    [lat, lon] = geographicGrid(raster, "gridvectors");

                    lat_in_range = lat(min(lats) <= lat & lat <= max(lats));
                    lon_in_range = lon(min(lons) <= lon & lon <= max(lons));
                    lon_lims = [min(lon_in_range, [], "all") max(lon_in_range, [], "all")];
                    lat_lims = [min(lat_in_range, [], "all") max(lat_in_range, [], "all")];

                    [~, raster] = geocrop(zeros(length(lat), length(lon)), raster, lat_lims, lon_lims);
                    row_idx = raster.LatitudeLimits(1) <= lat & lat <= raster.LatitudeLimits(2);
                    col_idx = raster.LongitudeLimits(1) <= lon & lon <= raster.LongitudeLimits(2);

                    indices = {row_idx, col_idx};
                    axes = {"lat", lat(row_idx), "lon", lon(col_idx)};
                otherwise
                    error("Unrecognized coordinate system type '%s'", raster_type);
            end
        end

        function [str] = format_filesize(sizes)
            % str = format_filesize(sizes)
            %   format byte sizes into <str> using SI prefixes up to TB

            arguments
                sizes {mustBeNumeric, mustBeNonempty};
            end

            pwr = floor(log10(sizes)/3)*3;
            units = dictionary([-Inf 0 3 6 9 12], ["B ", "B ", "KB", "MB", "GB", "TB"]);
            reduced_sizes = sizes./10.^pwr;
            reduced_sizes(~isfinite(reduced_sizes)) = 0;
            str = compose("%5.1f %s", reduced_sizes, units(pwr));
        end

        function [str] = sprintnf(str, name, rep)
            % str = sprintnf(str, Name, Value, ...)
            %   Adds name-based indexing to sprintf's %<number>$<spec>
            %   capability, so that %<name>$<spec> will use the corresponding
            %   Name-Value pair for formatting.
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
            str = sprintf(str, rep{:});
        end

        function [] = mustHaveCols(tbl, vars)
            % Internal utility
            % Validate that table has specified variables
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

        % function [] = mustHaveVars(mat, vars)
        %
        % end
    end
end

