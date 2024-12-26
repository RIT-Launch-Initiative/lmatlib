% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

% Methods of indexing:
% 

classdef (Abstract) nwpdata 
    properties (Constant, Access = private)
        metadata_cols = ["Element", "Unit", "Comment"];
    end
    
%% PUBLIC OBJECT METHODS
    methods (Static)
        function [data, raster, metadata] = read(path, params)
            arguments
                path (1,1) string {mustBeFile};
                params.elements (1,:) string;
                params.layers (1,:) string;
            end

            % Initialization
            elements = params.elements;
            layers = params.layers;
            info = georasterinfo(path);

            raster = info.RasterReference;
            raster_type = raster.CoordinateSystemType;
            switch raster_type
                case "planar"
                    [x, y] = worldGrid(raster, "gridvectors");
                    axes = {"y", y, "x", x};
                case "geographic"
                    [lat, lon] = geographicGrid(raster, "gridvectors");
                    axes = {"lat", lat, "lon", lon};
                otherwise
                    error("Unrecognized coordiante system type '%s'", raster_type);
            end

            metadata = info.Metadata;
            all_bands = nwpdata.find_bands(metadata, elements, layers);
            all_layers = categorical(metadata.ShortName(all_bands));
            all_layers = categories(all_layers, OutputType = "categorical");

            % keyboard;
            data_array = cell(1, length(elements));
            for i_element = 1:length(elements)
                element = elements(i_element);
                bands = nwpdata.find_bands(metadata, element, layers);
                [element_data, ~] = readgeoraster(path, Bands = bands);

                element_data = xarray(element_data, axes{:}, ...
                    layer = metadata.ShortName(bands), qty = element); 
                [~, reorder] = ismember(element_data.layer, all_layers);
                data_array{i_element} = element_data(:, :, reorder, :);
            end

            data = cat("qty", data_array{:});

            time = metadata.ValidTime(1);
            time.TimeZone = "UTC";
            data.time = time;

            if nargout == 3
                [~, i] = unique(metadata.Element);
                metadata = convertvars(metadata, "Element", "categorical");
                metadata = metadata(i, nwpdata.metadata_cols);
            end
        end
        
        function [data, raster, metadata] = read_baro(path, params)
            arguments
                path (1,1) string {mustBeFile};
                params.elements (1,:) string;
            end

            isbl_regex = "(\d+)-ISBL";
            [data, raster, metadata] = nwpdata.read(path, ...
                elements = params.elements, layers = isbl_regex);
            baro_levels = string(regexp(string(data.layer), isbl_regex, "tokens"));
            data.layer = str2double(baro_levels);
            data = sort(data, "layer", "descend");
        end

        function plot(data, raster, metadata)
            if any(size(data, ["time", "layer", "qty"]) > 1)
                error("Map data must be scalar along x/y or lat/lon dimensions");
            end
            ttl = sprintf("%s at layer %s\n%s", ...
                data.qty, data.layer, string(data.time, "yyyy/MM/dd HHmm z"));
            longname = metadata.Comment(metadata.Element == data.qty);

            raster_type = raster.CoordinateSystemType;
            switch raster_type
                case "planar"
                    ax = gca;
                    mapshow(ax, double(data), raster, DisplayType = "surface");
                case "geographic"
                    ax = axesm;
                    meshm(double(data), raster);
                otherwise
                    error("Unrecognized coordiante system type '%s'", raster_type);
            end
            title(ttl);
            cb = colorbar(ax);
            cb.Label.String = longname;
            axis(ax, "equal", "tight");
            
        end

        function filename = filename(model, date, forecast)
            arguments
                model (1,1) string {mustBeMember(model, "nam")};
                date (:, 1) datetime;
                forecast (:, 1) duration = hours(0);
            end

            % keyboard;
            % model_prefix = 
            filename = compose("nam_218_%s_%03d.grb2", string(date, "yyyyMMdd_HHmm"), ...
                floor(hours(forecast)));
        end
    end

%% INTERNAL UTILITIES
methods (Static, Access = private)
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

        % Get the URL for a NAM analysis

        %
        % function [filename, subfolder, url] = make_url(date)
        %     if any(mod(date.Hour, 6))
        %         error("Invalid model hour");
        %     end
        %
        %     base = "https://www.ncei.noaa.gov/data/north-american-mesoscale-model/access/analysis/";
        %     folder_fmt = "yyyyMM/yyyyMMdd/";
        %     file_fmt = "'nam_218'_yyyyMMdd_HHmm_000.'grb2'";
        %     subfolder = string(date, folder_fmt);
        %     filename = string(date, file_fmt);
        %     url = base + subfolder + filename;
        % end

        % function [name, url] = make_urls(product, date, forecast)
        %     arguments
        %         product (1,1) string;
        %         date (1,:) datetime;
        %         forecast (1, :) duration;
        %     end
        %     
        %     sources = dictionary;
        %     sources("nam") = struct(dir = "model-nam218")
        %     
        %     
        % end
        %
        % function [files, downloaded] = download(dest, product, dates, forecasts)
        %     [names, urls] = nwpdata.make_urls(product, dates, forecasts);
        %     files = fullfile(dest, names);
        %     
        %     downloaded = isfile(files);
        %     fprintf("Located %d files on disk:\n", sum(downloaded));
        %     fprintf(compose("\t%s\n", files(downloaded)));
        %
        %     on_server = isfile(urls(~downloaded));
        %     fprintf("Located %d files on server:\n", sum(on_server));
        %     
        %     if any(~present) 
        %         warning("Unable to locate %d files on server:\n%s", ...
        %             sum(~present), compose("\t%s\n", urls(~present)));
        %     end
        % end

        % function [files, downloaded] = download(dates, dest)
        %     [names, ~, urls] = nwpdata.make_url(dates);
        %     files = fullfile(dest, names);
        %     downloaded = false(size(files));
        %
        %     for i = 1:length(urls)
        %         fprintf("Downloading (%d) %s from %s\n", i, names(i), urls(i));
        %         if isfile(files(i))
        %             downloaded(i) = true;
        %             fprintf("Located on disk at %s\n", files(i));
        %         else
        %             tic;
        %             downloaded(i) = copyfile(urls(i), files(i));
        %             timed = toc;
        %             fprintf("Finished in %.2f sec\n", timed);
        %         end
        %
        %         if ~downloaded(i)
        %             warning("Unable to download from '%s'", names(i), urls(i));
        %         end
        %     end
        %     files = files(downloaded);
        % end
    end
end
