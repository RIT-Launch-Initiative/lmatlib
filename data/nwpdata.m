% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for reading and processing georeferenced rasters
% 
% Currently supported output products
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/gfs/">GFS (1.00, 0.50, 0.25-deg) pressure fields</a>
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/nam/">NAM (12, 3-km) pressure fields</a>
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/hrrr/">HRRR (3-km) pressure fields</a>
% Planned support
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/gens/">GEFS (0.50, 0.25-deg)</a>

classdef nwpdata < matlab.mixin.Scalar & handle & matlab.mixin.CustomDisplay
    properties (GetAccess = public, SetAccess = protected)
        model (1,1) string = missing;
        product (1,1) string = missing;
        cycle (1,1) datetime = NaT;
        forecast (1,:) duration;
        path (1,:) string;
        info (1,:) map.io.RasterInfo;
    end

    properties (Dependent)
        crs {mustBeA(crs, ["geocrs", "projcrs"])};
    end

    properties (Access = protected)
        model_id (1,1) string = missing;
        product_id (1,1) string = missing;
        filename (1,:) string;
        url (1,:) string;
    end
    
    properties (Constant, Access = protected)
        defs dictionary = nwpdata.data_definitions;
    end

    methods
        function obj = nwpdata(model, grid, cycle, forecasts)
            % ref = nwpdata(model, grid, cycle, forecasts)
            %   model   (string)    NWP model 
            %   grid    (string)    Data product
            %   cycle   (datetime)  Model output cycle (usually hourly or 6-hourly, 00 UTC)
            %   fcst    (duration)  Forecast extent (usually hourly or 3-hourly)
            arguments
                model (1,1) string;
                grid (1,1) string;
                cycle (1,1) datetime;
                forecasts (1,:) duration {mustBeNonempty};
            end 

            % Validate model
            if ~ismember(model, nwpdata.defs.keys)
                error("Model ""%s"" not recognized. Models supported: %s", ...
                    model, nwpdata.defs.keys.join(", "));
            end
            obj.model_id = model;
            model = nwpdata.defs{model};
            obj.model = model{"name"};

            % Validate product
            if ~ismember(grid, model{"products"}.keys)
                error("Grid ""%s"" not recognized. %s supports grids %s", ...
                    grid, upper(obj.model_id), model{"products"}.keys.join(", "));
            end
            product = model{"products"}{grid};
            obj.product_id = product{"wcoss_id"};
            obj.product = product{"name"};
            
            % Validate cycle
            if isempty(cycle.TimeZone)
                warning("nwpdata:timeZone", ...
                    "Model cycle does not have assigned time zone. Defaulting to UTC.");
            end
            cycle.TimeZone = "UTC";
            cycle_time = cycle - dateshift(cycle, "start", "day");
            if ~ismember(cycle_time, product{"cycle_values"})
                error("Invalid model cycle %s: %s-%s is produced %s", ...
                    cycle, model, grid, product{"cycle_hint"});
            end
            cycle.Format = "dd-MMM-yyyy HH z";
            obj.cycle = cycle;

            % Validate forecasts
            invalid = setdiff(forecasts, product{"forecast_values"});
            if ~isempty(invalid)
                error("Invalid model forecasts(s) %s: %s %s is produced %s", ...
                    mat2str(string(invalid)), upper(obj.model_id), grid, product{"forecast_hint"});
            end
            obj.forecast = forecasts;

            % Generate filenames and URLs
            fillin = @(fmt) nwpdata.populate(fmt, ...
                PRODUCT = obj.product_id, FF = hours(obj.forecast), FFF = hours(obj.forecast), ...
                YYYY = obj.cycle.Year, MM = obj.cycle.Month, DD = obj.cycle.Day, CC = obj.cycle.Hour);

            obj.filename = fillin(model{"filename"});
            src = fillin(model{"url"});
            obj.url = fullfile(src, obj.filename);
        end

        function download(obj, folder)
            arguments
                obj nwpdata;
                folder (1,1) string {mustBeFolder};
            end

            dest = fullfile(folder, string(obj.cycle, "yyyy-MM-dd"));
            [status, ~, ~] = mkdir(dest);
            if ~status
                error("Unable to create destination folder %s", dest);
            end

            [files, notfiles] = bulk_download(dest, obj.url, obj.filename);
            if ~isempty(notfiles)
                error("Unable to download %d file(s)", length(notfiles));
            end
            obj.path = files;

            for i_path = 1:length(obj.path)
                obj.info(i_path) = georasterinfo(obj.path(i_path));
            end
        end
        

        function crs = get.crs(obj)
            if isdownloaded(obj)
                crs = obj.info(1).CoordinateReferenceSystem;
            else 
                crs = [];
            end
        end

        function tru = isdownloaded(obj)
            tru = ~isempty(obj.path) && all(isfile(obj.path));
        end
        
        function data = read(obj, params)
            arguments
                obj nwpdata;   
                params.fields (1,:) string;
                params.layers (1,:) string;
                params.lats (1,2) double = [-Inf Inf];
                params.lons (1,2) double = [-Inf Inf];
            end

            downloadchk(obj);
            % NOTE: assumes all rasters are the same
            [~, indices, axes] = nwpdata.crop_raster(obj.info(end).RasterReference, ...
                params.lats, params.lons);
            fields = params.fields;
            layers = params.layers;

            metadata = obj.info(end).Metadata;
            all_bands = nwpdata.find_bands(metadata, fields, layers);
            output_fields = unique(metadata.Element(all_bands), "stable");
            output_layers = unique(metadata.ShortName(all_bands), "stable");
            times = obj.cycle + obj.forecast;

            % Allocate output
            nrows = length(axes{2});
            ncols = length(axes{4});
            ntimes = length(times);
            nlayers = length(output_layers);
            nfields = length(output_fields);

            data = xarray(NaN(nrows, ncols, ntimes, nlayers, nfields), ...
                axes{:}, time = times, layer = output_layers, field = output_fields);

            for i_time = 1:ntimes
                for i_layer = 1:nlayers
                    this_info = obj.info(i_time);
                    this_inventory = this_info.Metadata;
                    layer = output_layers(i_layer);
                    bands = nwpdata.find_bands(this_inventory, output_fields, layer);
                    fields = this_inventory.Element(bands);

                    if length(bands) < length(output_fields)
                        warning("nwpdata:missingFields", "Missing %s at layer %s", ...
                            mat2str(setdiff(output_fields, fields)), layer);
                    end

                    % [~, field_order] = ismember(fields, output_fields);
                    % Don't use <ismember>! 
                    % if only a subset of the fields exist, the data can't be assigned to using ':'
                    % This is the second time I've made this mistake
                    [present, i_output, i_file] = intersect(output_fields, fields, "stable");
                    if isempty(present)
                        continue;
                    end
                    layer_data = readgeoraster(this_info.Filename, Bands = bands);
                    data{:, :, i_time, i_layer, i_output} = ...
                        permute(layer_data(indices{:}, i_file), [1 2 4 5 3]);
                end
            end
        end
    end

    methods (Access = protected)
        function downloadchk(obj)
            err_id = "nwpdata:notDownloaded";
            if ~isdownloaded(obj)
                throwAsCaller(MException(err_id, "Files not downloaded." + ...
                    "Please run <data>.download(<folder>)"));
            end
        end

        function header = getHeader(~)
            header = "Referencing object for Numerical Weather Prediction data.";
            header = char(header);
        end

        function groups = getPropertyGroups(obj)
            propnames = string(properties(obj));
            proplist = struct;
            for i_prop = 1:length(propnames)
                propname = propnames(i_prop);
                prop = obj.(propname);
                if ~ismissing(prop) & ~isempty(prop)
                    proplist.(propname) = prop;
                end
            end
            
            groups = matlab.mixin.util.PropertyGroup(proplist);
        end
    end

%% INTERNAL UTILITIES
    methods (Static, Access = protected)
        function str = populate(str, name, rep)
            arguments
                str (:, 1) string;
            end
            arguments (Repeating)
                name (1, 1) string;
                rep (:, 1);
            end

            for i = 1:length(name)
                if isnumeric(rep{i})
                    fmt = sprintf("%%0%dd", strlength(name{i}));
                    rep{i} = compose(fmt, rep{i});
                end
                pat = "<" + name{i} + ">";
                if contains(str, pat)
                    str = strrep(str, pat, rep{i});
                end
            end
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

        function layers = find_bands(metadata, fields, layers)
            % [layers] = FIND_BANDS(metadata, fields, layers)
            % Find the layers (linear indices) matching any fields AND any layers
            %   The components of (fields) and (layers) are joined by
            %   "(e1|e2|...)" and treated as regular expressions
            arguments
                metadata table;
                fields (1, :) string;
                layers (1, :) string;
            end

            fields = "(" + join(fields, "|") + ")";
            layers = "(" + join(layers, "|") + ")";
            field_idx = ~cellfun(@isempty, regexp(metadata.Element, fields));
            layer_idx = ~cellfun(@isempty, regexp(metadata.ShortName, layers));
            layers = find(field_idx & layer_idx);
        end

        % Create weather data definitions
        function defs = data_definitions
            model_template = configureDictionary("string", "cell");
            model_template{"products"} = missing;
            model_template{"filename"} = missing;
            model_template{"url"} = missing;

            product_template = configureDictionary("string", "cell");
            product_template{"wcoss_id"} = missing;
            product_template{"name"} = missing;
            product_template{"grid_type"} = missing; % must be "geographic" or "planar"
            product_template{"cycle_values"} = [hours(0) hours(6) hours(12) hours(18)];
            product_template{"cycle_hint"} = "6-hourly at 00:00, 06:00, 12:00, 18:00";
            product_template{"forecast_values"} = NaT;
            product_template{"forecast_hint"} = NaT;

            % GFS definitions
            onedeg = product_template;
            onedeg{"wcoss_id"} = "pgrb2.1p00";
            onedeg{"name"} = "1.00 deg global latitude/longitude";
            onedeg{"grid_type"} = "geographic";
            onedeg{"forecast_values"} = hours(0:3:384);
            onedeg{"forecast_hint"} = "3-hourly up to 384 hours";

            halfdeg = onedeg;
            halfdeg{"wcoss_id"} = "pgrb2.0p50";
            halfdeg{"name"} = "0.50 deg global latitude/longitude";

            quartdeg = onedeg;
            quartdeg{"name"} = "0.25 deg global latitude/longitude";
            quartdeg{"forecast_values"} = hours(0:1:384);
            quartdeg{"forecast_hint"} = "hourly up to 384 hours";

            
            % USER FACING NAMES
            gfs_products = dictionary(["1.00 deg", "0.50 deg", "0.25 deg"], ...
                {onedeg, halfdeg, quartdeg});

            gfs = dictionary;
            gfs{"products"} = gfs_products;
            gfs{"name"} = "Global Forecast System";
            gfs{"filename"} = "gfs.t<CC>z.<PRODUCT>.f<FFF>";
            gfs{"url"} = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.<YYYY><MM><DD>/<CC>/atmos";

            % NAM definitions
            twelve = product_template;

            twelve{"wcoss_id"} = "awphys";
            twelve{"name"} = "12-km continental U.S. Lambert projected";
            twelve{"grid_type"} = "planar";
            twelve{"forecast_values"} = hours(0:1:84);
            twelve{"forecast_hint"} = "hourly up to 84 hours";

            three = product_template;
            three{"wcoss_id"} = "conusnest.hiresf";
            three{"name"} = "3-km contiguous U.S. Lambert projected grid";
            three{"grid_type"} = "planar";
            three{"forecast_values"} = hours(0:1:60);
            three{"forecast_hint"} = "hourly up to 60 hours";

            nam = model_template;
            nam{"products"} = dictionary(["12 km", "3 km"], {twelve three}); % USER FACING NAMES
            nam{"name"} = "North American Mesoscale";
            nam{"filename"} = "nam.t<CC>z.<PRODUCT><FF>.tm00.grib2";
            nam{"url"} = "https://noaa-nam-pds.s3.amazonaws.com/nam.<YYYY><MM><DD>";

            % HRRR definitions
            hrrr_product = product_template;
            hrrr_product{"wcoss_id"} = "wrfprsf";
            hrrr_product{"name"} = "3-km contiguous U.S. Lambert projected grid";
            hrrr_product{"grid_type"} = "planar";
            hrrr_product{"forecast_values"} = hours(0:1:60);
            hrrr_product{"forecast_hint"} = "hourly up to 48 hours";
            hrrr_product{"cycle_values"} = hours(0:1:23);
            hrrr_product{"cycle_hint"} = "hourly";

            hrrr = model_template;
            hrrr{"products"} = dictionary("3 km", {hrrr_product}); % USER FACING NAMES
            hrrr{"name"} = "High-Resolution Rapid Refresh";
            hrrr{"filename"} = "hrrr.t<CC>z.<PRODUCT><FF>.grib2";
            hrrr{"url"} = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.<YYYY><MM><DD>/conus";

            % USER FACING NAMES
            defs = dictionary(["gfs", "nam", "hrrr"], {gfs, nam, hrrr}); 
        end
    end
end
