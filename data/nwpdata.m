% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

classdef (Sealed) nwpdata
    properties (GetAccess = public, SetAccess = private)
        data (1, :, :, :, :) double; % t, y, x, p, d:
        raster map.rasterref.MapPostingsReference;
        isobars (1, :) double;
        metadata table;
    end

    properties (Dependent)
        % x
        % y
        % z
    end

    methods
        function obj = nwpdata(path, elements, latlon, side, levels)
            arguments
                path (1,1) string;
                elements (1,:) string;
                latlon (1,2) double = [0, 0];
                side (1,1) double = Inf;
                levels (1,1) string = "";
            end
            
            info = georasterinfo(path);
            metadata = info.Metadata;
            raster = info.RasterReference;
            [obj.raster, x_idx, y_idx] = nwpdata.cropraster(raster, latlon(1), latlon(2), side);

            if levels == "" % all isobaric levels
                levels = "^(\d+)(?=-ISBL$)";
                ptn = regexpPattern(levels);

                % Find all unique isobaric levels 
                all_bands = nwpdata.find_bands(metadata, elements, levels);
                isobars = extract(metadata{all_bands, "ShortName"}, ptn);
                isobars = str2double(isobars);
                obj.isobars = sort(unique(isobars), "descend");

                % Now enough is known to allocate the object data member
                obj.data = NaN([1, obj.raster.RasterSize, length(obj.isobars), length(elements)]);

                for i_element = 1:length(elements)
                    % Find all isobaric levels for current data element
                    bands = nwpdata.find_bands(metadata, elements(i_element), levels);
                    isobars = str2double(extract(metadata{bands, "ShortName"}, ptn));

                    [data, ~] = readgeoraster(path, Bands = bands);

                    % "Shouldn't" be necessary to intersect and match indices
                    % since the data appears to be consistently sorted, but
                    % just to be safe
                    [~, obj_idx, data_idx] = intersect(obj.isobars, isobars, "stable");
                    % Put data into object
                    obj.data(1, :, :, obj_idx, i_element) = data(y_idx, x_idx, data_idx);
                end
            else
                bands = nwpdata.find_bands(metadata, elements, levels);
                obj.data = NaN([1, obj.raster.RasterSize, 1, length(bands)]);

                [data, ~] = readgeoraster(path, Bands = bands);

                obj.data(1, :, :, 1, :) = data(y_idx, x_idx, :);
                obj.metadata = metadata;
                obj.isobars = [];
            end
        end

        % function tab = latlon(obj, lat, lon, method)
        %     arguments
        %         obj nwpdata;
        %         lat (1,1) double;
        %         lon (1,1) double;
        %         method (1,1) string {check_strings(method, ["interpolate", "closest"])} = "interpolate";
        %     end
        % end
    end

    methods (Static)
        function [raster, x_idx, y_idx] = cropraster(raster, lat, lon, side)
            if side == Inf
                x_idx = true(1, raster.RasterSize(2));
                y_idx = true(1, raster.RasterSize(1));
                return;
            end

            [x_values, y_values] = worldGrid(raster, "gridvectors");
            [x_center, y_center] = projfwd(raster.ProjectedCRS, lat, lon);
            x_limits = x_center + side/2*[-1 1];
            y_limits = y_center + side/2*[-1 1];
            [~, raster] = mapcrop(zeros(raster.RasterSize), raster, x_limits, y_limits);
            x_idx = raster.XWorldLimits(1) <= x_values & x_values <= raster.XWorldLimits(2);
            y_idx = raster.YWorldLimits(1) <= y_values & y_values <= raster.YWorldLimits(2);

            assert(sum(x_idx) == raster.RasterSize(2), "X indexor must match raster size");
            assert(sum(y_idx) == raster.RasterSize(1), "Y indexor must match raster size");
        end

        function bands = find_bands(metadata, elements, levels)
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
