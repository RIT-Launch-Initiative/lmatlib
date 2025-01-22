% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for reading and processing georeferenced rasters
% 
% Currently supported output products
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/gfs/">GFS (1.00, 0.50, 0.25-deg) pressure fields</a>
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/nam/">NAM (12, 3-km) pressure fields</a>
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/hrrr/">HRRR (3-km) pressure fields</a>
% Planned support
% - <a href="https://www.nco.ncep.noaa.gov/pmb/products/gens/">GEFS (0.50, 0.25-deg)</a>

classdef nwpdata 


    properties (GetAccess = public, SetAccess = protected)
        % member (1,1) {mustBeA(member, ["string", "numeric"])};
        model_name (1,1) string;
        product_name (1,1) string;
        file (1,1) string = missing;
    end

    properties (Dependent)
        time (1,1) datetime;
        inventory table;
    end

    properties (Hidden)
        model_id (1,1) string = missing;
        product (1,1) string = missing;
        cycle (1,1) datetime = missing;
        forecast (1,:) duration = missing;
        url (1,1) string = missing;
        info (1,1) map.io.RasterInfo = missing;
    end
    
    properties (Constant, Access = protected)
        defs dictionary = nwpdata.data_definitions;
    end

    methods
        function obj = nwpdata(model, grid, cycle, forecasts)
            arguments
                model (1,1) string;
                grid (1,1) string;
                cycle (1,1) datetime;
                forecasts (1,:) duration;
            end 

            mustBeMember(model, nwpdata.defs.keys);
            obj.model_id = model;
            obj.model_name = nwpdata.defs{model}{"name"};

            products = nwpdata.defs{model}{"products"};
            mustBeMember(grid, products.keys);
            obj.product = products{grid}{"product"};
            obj.product_name = products{grid}{"name"};
            
            cycle.TimeZone = "UTC";
            cycle_time = cycle - dateshift(cycle, "start", "day");
            if ~ismember(cycle_time, products{grid}{"cycle_values"})
                error("Invalid model cycle %s: %s-%s is produced %s", ...
                    cycle, model, grid, products{grid}{"cycle_hint"});
            end
            obj.cycle = cycle;

            invalid = setdiff(forecasts, products{grid}{"forecast_values"});
            if ~isempty(invalid)
                error("Invalid model forecasts(s) %s: %s-%s is produced %s", ...
                    mat2str(string(invalid)), model, grid, products{grid}{"forecast_hint"});
            end
            obj.forecast = forecasts;
        end

        % function files = download(obj, folder)
        %     
        % end
        %
        % function output = read(obj, params)
        %     
        % end

        function time = get.time(obj)
            time = obj.cycle + obj.forecast;
            % fcst_text = string(obj.cycle, sprintf("(dd-MMM-yyyy HH z'+%g')", hours(obj.forecast)));
            % time.Format = "dd-MMM-yyyy HHz" + "'" + fcst_text + "'";
        end

        function inv = get.inventory(obj)
            if ismissing(obj.info)
                inv = missing;
            else
                inv = obj.info.Metadata;
            end
        end
    end
    methods (Static)

        function defs = data_definitions
            model_template = configureDictionary("string", "cell");
            model_template{"products"} = missing;
            model_template{"filename"} = missing;
            model_template{"url"} = missing;


            product_template = configureDictionary("string", "cell");
            product_template{"product"} = missing;
            product_template{"name"} = missing;
            product_template{"grid_type"} = missing; % must be "geographic" or "planar"
            product_template{"cycle_values"} = [hours(0) hours(6) hours(12) hours(18)];
            product_template{"cycle_hint"} = "6-hourly at 00:00, 06:00, 12:00, 18:00";
            product_template{"forecast_values"} = NaT;
            product_template{"forecast_hint"} = NaT;
            % product_template{"member_values"} = NaN;
            % product_template{"member_hint"} = "Not an ensemble forecast";

            % GFS definitions
            onedeg = product_template;
            onedeg{"product"} = "pgrb2.1p00";
            onedeg{"name"} = "1.00 deg global latitude/longitude";
            onedeg{"grid_type"} = "geographic";
            onedeg{"forecast_values"} = hours(0:3:384);
            onedeg{"forecast_hint"} = "3-hourly up to 384 hours";

            halfdeg = onedeg;
            halfdeg{"product"} = "pgrb2.0p50";
            halfdeg{"name"} = "0.50 deg global latitude/longitude";

            quartdeg = onedeg;
            quartdeg{"name"} = "0.25 deg global latitude/longitude";
            quartdeg{"forecast_values"} = hours(0:1:384);
            quartdeg{"forecast_hint"} = "hourly up to 384 hours";

            gfs_products = dictionary(["1.00 deg", "0.50 deg", "0.25 deg"], {onedeg, halfdeg, quartdeg});

            gfs = dictionary;
            gfs{"products"} = gfs_products;
            gfs{"name"} = "Global Forecast System";
            gfs{"filename"} = "gfs.t<CC>z.<PRODUCT>.f<FFF>";
            gfs{"url"} = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.<YYYY><MM><DD>/<CC>/atmos";

            % NAM definitions
            twelve = product_template;

            twelve{"product"} = "awphys";
            twelve{"name"} = " 12-km continental U.S. Lambert projected";
            twelve{"grid_type"} = "planar";
            twelve{"forecast_values"} = hours(0:1:84);
            twelve{"forecast_hint"} = "hourly up to 84 hours";

            three = product_template;
            three{"product"} = "conusnest.hiresf";
            three{"name"} = "3-km contiguous U.S. Lambert projected grid";
            three{"grid_type"} = "planar";
            three{"forecast_values"} = hours(0:1:60);
            three{"forecast_hint"} = "hourly up to 60 hours";

            nam = model_template;
            nam{"products"} = dictionary(["12 km", "3 km"], {twelve three});
            nam{"name"} = "North American Mesoscale";
            nam{"filename"} = "nam.t<CC>z.<PRODUCT><FF>.tm00.grib2";
            nam{"url"} = "https://noaa-nam-pds.s3.amazonaws.com/nam.<YYYY><MM><DD>";

            % HRRR definitions
            hrrr_product = product_template;
            hrrr_product{"product"} = "conusnest.hiresf";
            hrrr_product{"name"} = "3-km contiguous U.S. Lambert projected grid";
            hrrr_product{"grid_type"} = "planar";
            hrrr_product{"forecast_values"} = hours(0:1:60);
            hrrr_product{"forecast_hint"} = "hourly up to 48 hours";
            hrrr_product{"cycle_values"} = hours(0:1:23);
            hrrr_product{"cycle_hint"} = "hourly";

            hrrr = model_template;
            hrrr{"product"} = dictionary("3 km", hrrr_product);
            hrrr{"name"} = "High-Resolution Rapid Refresh";
            hrrr{"filename"} = "hrrr.t<CC>z.<PRODUCT><FF>.grib2";
            hrrr{"url"} = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.<YYYY><MM><DD>/conus";

            defs = dictionary(["gfs", "nam", "hrrr"], {gfs, nam, hrrr}); 
        end
    end

%% PUBLIC OBJECT METHODS
    methods (Static)
        function [data, raster, metadata] = read(path, params)
            % Read arbitrary data from GRIB2 file into xarray
            % [data, raster, metadata] = nwpdata.read(path, Name = Value)
            %   INPUTS
            %   path        Path to GRIB2 file
            %   
            %   Required name-value arguments
            %   fields      List of data fields (TMP, HGT, ...)
            %               May be regex
            %   layers      List of data layers (0-SFC, 100000-ISBL, EATM)
            %               May be regex
            %   
            %   Optional name-value arguments
            %   lats        Latitude limits 
            %               default [-Inf Inf]
            %   lons        Longitude limits 
            %               default [-Inf Inf]
            %   
            %   OUTPUTS
            %   data        xarray with dimensions (y/x or lat/lon, layer, field)
            %   raster      Location referencing object
            %   metadata    Table containing information about layers and fields
            arguments
                path (1,1) {mustBeFile}
                params.fields (1,:) string;
                params.layers (1,:) string;

                params.lats (1,2) double = [-Inf Inf];
                params.lons (1,2) double = [-Inf Inf];
            end

            % Initialization
            fields = params.fields;
            layers = params.layers;
            info = georasterinfo(path);

            raster = info.RasterReference;
            [raster, indices, axes] = nwpdata.crop_raster(raster, params.lats, params.lons);

            metadata = info.Metadata;
            all_bands = nwpdata.find_bands(metadata, fields, layers);
            output_fields = unique(metadata.Element(all_bands));
            output_layers = unique(metadata.ShortName(all_bands));

            nrows = length(axes{2});
            ncols = length(axes{4});
            nfields = length(output_fields);
            nlayers = length(output_layers);

            % Allocate output
            data = xarray(NaN(nrows, ncols, nlayers, nfields), ...
                axes{:}, layer = output_layers, field = output_fields);

            for i_layer = 1:nlayers
                layer = output_layers(i_layer);
                bands = nwpdata.find_bands(metadata, output_fields, layer);
                fields = metadata.Element(bands);

                if length(bands) < length(output_fields)
                    warning("Missing %s at layer %s", ...
                        layer, mat2str(setdiff(output_fields, fields)));
                end

                [~, field_order] = ismember(fields, output_fields);
                layer_data = readgeoraster(info.Filename, Bands = bands);
                data{:, :, i_layer, field_order} = permute(layer_data(indices{:}, :), [1 2 4 3]);
            end

            time = metadata.ValidTime(1);
            time.TimeZone = "UTC";
            data.time = time;

            if nargout == 3
                metadata = metadata(all_bands);
            end
        end

        % function data = read_internal(info, fields, layers)
        %     arguments
        %         info map.io.RasterInfo
        %         fields (1, :) string;
        %         layers (1, :) string;
        %     end 
        %
        %     metadata = info.Metadata;
        %     bands = nwpdata.find_bands(metadata, fields, layers);
        %     layers = unique(metadata.ShortName(bands));
        % end

        function [filename, url] = filename(model, product, date, forecast)
            arguments
                model (1,1) string {mustBeMember(model, ["nam", "gfs", "hrrr"])};
                product (1,1) string;
                date (:, 1) datetime;
                forecast (:, 1) duration = hours(0);
            end 
            
            if ~isscalar(date) && ~isscalar(forecast) ...
                && length(date) ~= length(forecast)
                error("Date and forecast must be the same length if both are nonscalar");
            end

            % Define format replacements
            reps = dictionary;
            reps{"<FF>"} = compose("%02d", hours(forecast));
            reps{"<FFF>"} = compose("%03d", hours(forecast));
            reps{"<YYYY>"} = compose("%04d", date.Year);
            reps{"<MM>"} = compose("%02d", date.Month);
            reps{"<DD>"} = compose("%02d", date.Day);
            reps{"<CC>"} = compose("%02d", date.Hour);
            reps{"<PRODUCT>"} = product;

            % https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.20140731/conus/hrrr.t00z.wrfnatf01.grib2
            switch model
                case "gfs"
                    filename = "gfs.t<CC>z.<PRODUCT>.f<FFF>";
                    folder = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.<YYYY><MM><DD>/<CC>/atmos";
                case "nam"
                    filename = "nam.t<CC>z.<PRODUCT><FF>.tm00.grib2";
                    folder = "https://noaa-nam-pds.s3.amazonaws.com/nam.<YYYY><MM><DD>";
                case "hrrr"
                    filename = "hrrr.t<CC>z.<PRODUCT><FF>.grib2";
                    folder = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.<YYYY><MM><DD>/conus";
            end

            keys = reps.keys;
            for i = 1:length(keys)
                key = keys(i);
                value = reps{key};
                filename = strrep(filename, key, value);
                folder = strrep(folder, key, value);
            end
            url = fullfile(folder, filename);
        end

        function str = populate(str, name, rep)
            arguments
                str (:, 1) string;
            end
            arguments (Repeating)
                name (1, 1) string;
                rep (:, 1) string;
            end

            for i = 1:length(name)
                if isnumeric(rep{i})
                    fmt = sprintf("%%0%dd", strlength(name{i}));
                    rep{i} = compose(fmt, rep{i});
                end
                str = strrep(str, "<" + name{i} + ">", rep{i});
            end
        end

        function files = download(dest, model, product, date, forecast)
            arguments
                dest (1,1) string {mustBeFolder};
                model (1,1) string;
                product (1,1) string;
                date (:, 1) datetime;
                forecast (:, 1) duration;
            end 
            
            % create paths
            dest = fullfile(dest, string(date, "yyyy-MM-dd"));
            [status, ~, ~] = mkdir(dest);
            if ~status
                error("Unable to create %s", dest);
            end

            [filename, url] = nwpdata.filename(model, product, date, forecast);
            files = bulk_download(dest, url, filename);
        end

        function coords = latlon2raster(raster, lat, lon)
            % Project geographic coordiantes onto the raster's coordinate system
            arguments
                raster
                lat double;
                lon double;
            end

            type = raster.CoordinateSystemType;
            
            switch type
                case "planar"
                    [coords.x, coords.y] = projfwd(raster.ProjectedCRS, lat, lon);
                case "geographic"
                    coords.lat = lat;
                    coords.lon = lon;
                otherwise
                    error("Unrecognized coordinate system type '%s'", raster_type);
            end
        end
    end

%% INTERNAL UTILITIES
methods (Static, Access = protected)
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
    end
end
