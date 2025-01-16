% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

% Methods of indexing:
% 

classdef (Abstract) nwpdata 
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
            %   layers      List of data layers (0-SFC, 100000-ISBL, EATM)
            %   
            %   Optional name-value arguments are used to crop the read data to
            %   a square to save memory. For planar projected maps, lat and lon
            %   are projected onto the map to get y/x limits. For geographic
            %   maps, side() is converted to angle limits using the Earth's
            %   radius, and the mean of the latitude limits is used to find the
            %   angle limits in longitude
            %   lat         Latitude enter of crop
            %   lon         Longitude of center of crop
            %   side        Side length of square used to crop
            %   
            %   OUTPUTS
            %   data        xarray with dimensions (y/x or lat/lon, layer, field)
            %   raster      Location referencing object
            %   metadata    Table containing information about layers and fields
            
            arguments
                path (1,1) string {mustBeFile};
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
            all_layers = metadata.ShortName(all_bands);
            layer_template = unique(all_layers);
            data_array = cell(1, length(fields));

            for i_field = 1:length(fields)
                field = fields(i_field);
                bands = nwpdata.find_bands(metadata, field, layers);
                layers = metadata.ShortName(bands);

                if length(layers) ~= length(layer_template)
                    error("Missing data for field %s at %s", field, ...
                        mat2str(setdiff(layer_template, layers)));
                end

                [field_data, ~] = readgeoraster(path, Bands = bands);
                [~, ~, order] = intersect(layer_template, layers, "stable");
                data_array{i_field} = field_data(indices{:}, order);
            end

            args = {};
            if ~isscalar(layer_template)
                args = [args {"layer", layer_template}];
            end
            if ~isscalar(fields)
                args = [args {"field", fields}];
            end
            data = xarray(cat(4, data_array{:}), ...
                axes{:}, args{:});

            if nargout == 3
                metadata = metadata(all_bands, ["ShortName", "Element", "Unit", "Comment", "ValidTime"]);
                metadata = convertvars(metadata, ["ShortName", "Element"], "categorical");
            end
        end
        
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

        function [x_out, y_out] = outline(x, y)
            arguments
                x (1, :) double;
                y (1, :) double;
            end
            x_out = [repmat(x(1), 1, length(y)), x, repmat(x(end), 1, length(y)), flip(x)];
            y_out = [y, repmat(y(end), 1, length(x)), flip(y), repmat(y(1), 1, length(x))];
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
