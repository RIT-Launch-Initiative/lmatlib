% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

classdef (Sealed) nwpdata
    properties (GetAccess = public, SetAccess = private)
        data (1, :, :, :, :) double; % t, y, x, p, d:
        time (1,1) datetime;
        raster map.rasterref.MapPostingsReference;
        isobars (1, :) double;
        metadata table;

        created (1,1) datetime;
    end

    % indexing - go into data() and cut time, raster, metadata, isobars accordingly

    % properties (Dependent)
    %     grid_vectors
    %     x_limits
    %     y_limits
    %     xy_grids
    % end

    methods
        function obj = nwpdata(path, params)
            arguments
                path (1,1) string;
                params.elements (1,:) string;
                params.lat (1,1) double = 0;
                params.lon (1,1) double = 0;
                params.side (1,1) double = Inf;
                params.levels (1,:) = [-Inf Inf];
            end

            % Initialization
            elements = params.elements;
            levels = params.levels;
            
            info = georasterinfo(path);
            metadata = info.Metadata;
            raster = info.RasterReference;
            [obj.raster, x_idx, y_idx] = nwpdata.cropraster(raster, ...
                params.lat, params.lon, params.side);

            if isnumeric(levels)
                ptn = "^(\d+)-ISBL";
                all_bands = nwpdata.find_bands(metadata, elements, ptn);
                isobars = string(regexp(metadata.ShortName(all_bands), ptn, "tokens"));
                isobars = str2double(isobars);
                % cut levels to specified range
                isobars = isobars(min(levels) <= isobars & isobars <= max(levels));

                obj.isobars = sort(unique(isobars), "descend");

                % Now enough is known to allocate the object data member
                obj.data = NaN([1, obj.raster.RasterSize, ...
                    length(obj.isobars), length(elements)]);

                % read elements in one at a time - all-at-once takes too much memory
                for i_element = 1:length(elements)
                    % Find all isobaric levels for current data element
                    bands = nwpdata.find_bands(metadata, elements(i_element), ptn);
                    isobars = string(regexp(metadata.ShortName(bands), ptn, "tokens"));
                    isobars = str2double(isobars);
                    % Not a fan of repeating this code, but don't have a better idea
                    in_range = min(levels) <= isobars & isobars <= max(levels);
                    isobars = isobars(in_range);
                    bands = bands(in_range);

                    [data, ~] = readgeoraster(path, Bands = bands);

                    % "Shouldn't" be necessary to intersect and match indices
                    % since the data appears to be consistently sorted, but
                    % just to be safe
                    [~, obj_idx, data_idx] = intersect(obj.isobars, isobars, "stable");
                    if length(obj_idx) < size(obj.data, 4)
                        warning("Missing data for element %s at %d isobars", ...
                            elements(i_element), size(obj.data, 4) - length(obj_idx));
                    end
                    % Put data into object
                    obj.data(1, :, :, obj_idx, i_element) = data(y_idx, x_idx, data_idx);
                    obj.metadata = [obj.metadata; metadata(bands(1), :)];
                end
            else
                % Read the bands and report them
                bands = nwpdata.find_bands(metadata, elements, levels);
                obj.data = NaN([1, obj.raster.RasterSize, 1, length(bands)]);

                [data, ~] = readgeoraster(path, Bands = bands);

                obj.data(1, :, :, 1, :) = data(y_idx, x_idx, :);
                obj.metadata = metadata(bands, :);
                obj.isobars = [];
            end
            obj.time = metadata.ValidTime(1);
            obj.created = metadata.ReferenceTime(1);
            obj.metadata = obj.metadata(:, ["Element", "Unit", "Comment"]);
        end
    end

    methods (Static)
        function [raster, x_idx, y_idx] = cropraster(raster, lat, lon, side)
            % [raster, x_idx, y_idx] = CROPRASTER(raster, lat, lon, side)
            % Crop the (raster) object to a square centered at (lat, lon) that is (side) wide
            % Return the logical index vectors to crop the associated matrix
            % (separated out because they are used multiple times)
            if side == Inf
                x_idx = true(1, raster.RasterSize(2));
                y_idx = true(1, raster.RasterSize(1));
                return;
            end

            [x_values, y_values] = worldGrid(raster, "gridvectors");
            [x_center, y_center] = projfwd(raster.ProjectedCRS, lat, lon);
            x_limits = x_center + side/2*[-1 1];
            y_limits = y_center + side/2*[-1 1];
            if isempty(raster)
                error("Cropping produced empty result");
            end
            [~, raster] = mapcrop(zeros(raster.RasterSize), raster, x_limits, y_limits);
            x_idx = raster.XWorldLimits(1) <= x_values & x_values <= raster.XWorldLimits(2);
            y_idx = raster.YWorldLimits(1) <= y_values & y_values <= raster.YWorldLimits(2);

            assert(sum(x_idx) == raster.RasterSize(2), "X indexor must match raster size");
            assert(sum(y_idx) == raster.RasterSize(1), "Y indexor must match raster size");
        end

        function bands = find_bands(metadata, elements, levels)
            % [bands] = FIND_BANDS(metadata, elements, levels)
            % Find the bands (linear indices) matching any elements AND any levels
            %   The components of (elements) and (levels) are joined by
            %   "(e1|e2|...)" and treated as regular expressions
            arguments
                metadata table;
                elements (1, :) string;
                levels (1, :) string;
            end
            elements = "(" + join(elements, "|") + ")";
            levels = "(" + join(levels, "|") + ")";
            element_idx = ~cellfun(@isempty, regexp(metadata.Element, elements));
            level_idx = ~cellfun(@isempty, regexp(metadata.ShortName, levels));
            bands = find(element_idx & level_idx);
        end

        % Get the URL for a NAM analysis
        function [filename, subfolder, url] = make_url(date)
            % https://www.ncei.noaa.gov/data/north-american-mesoscale-model/access/analysis/
            % https://www.ncei.noaa.gov/data/north-american-mesoscale-model/access/analysis/202112/20211222/nam_218_20211222_0000_001.grb2
            base = "https://www.ncei.noaa.gov/data/north-american-mesoscale-model/access/analysis/";
            folder_fmt = "yyyyMM/yyyyMMdd/";
            file_fmt = "'nam_218'_yyyyMMdd_HHmm_000.'grb2'";
            subfolder = string(date, folder_fmt);
            filename = string(date, file_fmt);
            url = base + subfolder + filename;
        end
    end

    methods (Static, Access = private)
        function check_strings(input, strings, prefix)
            % Validate that all members of the input are in the defined strings
            arguments
                input string;
                strings string;
                prefix (1,1) string = "";
            end
            isn_t = find(~ismember(input, strings));
            if isempty(isn_t)
                return;
            end

            string_print = compose("\t%s\n", strings).join("");
            if isscalar(input)
                text = "Input must be one of";
            elseif isscalar(isn_t)
                text = sprintf("Input at index %d is '%s', must be one of", input(isn_t), isn_t);
            else
                text = sprintf("Input at index %s must be one of", mat2str(isn_t));
            end

            ME = MException("openrocket:invalid_value", sprintf("%s%s\n%s", prefix, text, string_print));
            throwAsCaller(ME);
        end
    end
end
