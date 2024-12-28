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
            %   elements    List of data elements (TMP, HGT, ...)
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
            %   data        xarray with dimensions (y/x or lat/lon, layer, element, time)
            %   raster      Location referencing object
            %   metadata    Table containing information about layers and elements
            
            arguments
                path (1,1) string {mustBeFile};
                params.elements (1,:) string;
                params.layers (1,:) string;

                params.lat (1, 1) double = 0;
                params.lon (1, 1) double = 0;
                params.side (1, 1) double = Inf;
            end

            % Initialization
            elements = params.elements;
            layers = params.layers;
            info = georasterinfo(path);

            raster = info.RasterReference;
            [raster, indices, axes] = nwpdata.crop_raster(raster, params.lat, params.lon, params.side);

            metadata = info.Metadata;
            all_bands = nwpdata.find_bands(metadata, elements, layers);
            all_layers = metadata.ShortName(all_bands);
            layer_template = unique(all_layers);
            data_array = cell(1, length(elements));

            for i_element = 1:length(elements)
                element = elements(i_element);
                bands = nwpdata.find_bands(metadata, element, layers);
                layers = metadata.ShortName(bands);
                if length(layers) ~= length(layer_template)
                    error("Missing data for element %s at %s", element, ...
                        mat2str(setdiff(layer_template, layers)));
                end

                [element_data, ~] = readgeoraster(path, Bands = bands);
                [~, ~, order] = intersect(layer_template, layers, "stable");
                data_array{i_element} = element_data(indices{:}, order);
            end

            data = xarray(cat(4, data_array{:}), ...
                axes{:}, layer = layer_template, element = elements);

            time = metadata.ValidTime(1);
            time.TimeZone = "UTC";
            data.time = time;

            if nargout == 3
                metadata = metadata(all_bands, ["ShortName", "Element", "Unit", "Comment"]);
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
            files = fullfile(dest, filename);
            
            % check presence of every file
            on_disk = false(size(files));
            on_server = false(size(files));
            sizes = zeros(size(files));
            for i_file = 1:length(files)
                info = dir(files(i_file));
                if isscalar(info)
                    on_disk(i_file) = true;
                    sizes(i_file) = info.bytes;
                    continue;
                end
                
                info = dir(url(i_file));
                if isscalar(info)
                    on_server(i_file) = true;
                    sizes(i_file) = info.bytes;
                    continue;
                end
            end
            
            % emit status messages for files located on disk, on server, or not at all
            n_on_disk = sum(on_disk);
            sz_on_disk = sum(sizes(on_disk));
            files_on_disk = files(on_disk);
            if n_on_disk > 0
                fprintf("Located %d files on disk (%s):\n%s", n_on_disk, format_file_size(sz_on_disk), ...
                    long_list(files_on_disk));
            end

            n_on_server = sum(on_server);
            sz_on_server = sum(sizes(on_server));
            files_on_server = url(on_server);
            if n_on_server > 0
                fprintf("Located %d files on server (%s):\n%s", n_on_server, format_file_size(sz_on_server), ...
                    long_list(files_on_server));
            end

            n_nowhere = sum(~on_disk & ~on_server);
            files_nowhere = url(~on_disk & ~on_server);
            if n_nowhere > 0
                warning("%d files not found: \n%s", n_nowhere, long_list(files_nowhere));
            end

            % if sz_on_server > 500e6
            %     prompt = sprintf("Download %s [y/n]: ", format_file_size(sz_on_server));
            %     res = input(prompt, "s");
            %     if res ~= "y"
            %         files = files(on_disk);
            %         return;
            %     end
            % end

            % download the files located on server
            for i_file = find(on_server(:)')
                file_url = url(i_file);
                fprintf("Downloading %s from %s... ", format_file_size(sizes(i_file)), file_url)
                tic;
                downloaded = websave(files(i_file), file_url);
                % downloaded = copyfile(file_url, dest);
                timed = toc;
                if downloaded
                    on_disk(i_file) = true;
                    fprintf("succeeded in %.1f sec\n", timed);
                else
                    fprintf("failed\n")
                end
            end

            files = files(on_disk);
            
            % Helper functions
            function str = format_file_size(sz)
                pwr = floor(log10(sz)/3)*3;
                units = ["B", "KB", "MB", "GB", "TB"];
                logs = [0 3 6 9 12];


                prefix = units(logs == pwr);
                str = sprintf("%0.2f %s", sz / 10^pwr, prefix);
            end

            function str = long_list(strs)
                str = compose("\t%s\n", strs).join("");
            end
        end
    end

%% INTERNAL UTILITIES
methods (Static, Access = private)
        function [raster, indices, axes] = crop_raster(raster, lat_c, lon_c, side)
            raster_type = raster.CoordinateSystemType;
            switch raster_type
                case "planar"
                    [x, y] = worldGrid(raster, "gridvectors");

                    % if we need to crop, convert parameters to raster's coordinate system
                    if isfinite(side)
                        [x_c, y_c] = projfwd(raster.ProjectedCRS, lat_c, lon_c);
                        x_lims = x_c + side/2*[-1 1];
                        y_lims = y_c + side/2*[-1 1];

                        [~, raster] = mapcrop(zeros(length(y), length(x)), raster, x_lims, y_lims);
                        % need to use the new limits instead of x_lims/y_lims because mapcrop
                        % leaves an extra element on either side
                        row_idx = raster.YWorldLimits(1) <= y & y <= raster.YWorldLimits(2);
                        col_idx = raster.XWorldLimits(1) <= x & x <= raster.XWorldLimits(2);
                        
                        indices = {row_idx, col_idx};
                        axes = {"y", y(row_idx), "x", x(col_idx)};
                    else % pick everything without doing math
                        indices = {':', ':'};
                        axes = {"y", y, "x", x};
                    end
                case "geographic"
                    [lat, lon] = geographicGrid(raster, "gridvectors");
                    
                    % if we need to crop, convert parameters to raster's coordinate system
                    if isfinite(side)
                        R = raster.GeographicCRS.Spheroid.MeanRadius;
                        % convert length to degrees
                        lat_lims = lat_c + rad2deg(side/R)/2 * [-1 1];
                        lon_lims = lon_c + rad2deg(side/(R*cosd(mean(lat_lims))))/2 * [-1 1];
                        
                        [~, raster] = geocrop(zeros(length(lat), length(lon)), raster, lat_lims, lon_lims);
                        % use new limits, same reason as for mapcrop() above
                        row_idx = raster.LatitudeLimits(1) <= lat & lat <= raster.LatitudeLimits(2);
                        col_idx = raster.LongitudeLimits(1) <= lon & lon <= raster.LongitudeLimits(2);

                        indices = {row_idx, col_idx};
                        axes = {"lat", lat(row_idx), "lon", lon(col_idx)};
                    else % pick everything without doing mat
                        indices = {':', ':'};
                        axes = {"lat", lat, "lon", lon};
                    end
                otherwise
                    error("Unrecognized coordiante system type '%s'", raster_type);
            end
        end

        function layers = find_bands(metadata, elements, layers)
            % [layers] = FIND_BANDS(metadata, elements, layers)
            % Find the layers (linear indices) matching any elements AND any layers
            %   The components of (elements) and (layers) are joined by
            %   "(e1|e2|...)" and treated as regular expressions
            arguments
                metadata table;
                elements (1, :) string;
                layers (1, :) string;
            end

            elements = "(" + join(elements, "|") + ")";
            layers = "(" + join(layers, "|") + ")";
            element_idx = ~cellfun(@isempty, regexp(metadata.Element, elements));
            layer_idx = ~cellfun(@isempty, regexp(metadata.ShortName, layers));
            layers = find(element_idx & layer_idx);
        end
    end
end
