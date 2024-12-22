% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

classdef (Sealed) nwpdata < matlab.mixin.indexing.RedefinesParen

    properties (GetAccess = public, SetAccess = private)
        data (:, :, :, :, :) double; 
        time (1,:) datetime;
        levels (1, :);
        raster map.rasterref.MapPostingsReference;
        metadata table;
        created (1,:) datetime;
    end

    properties (Constant, Access = private)
        press_regex = "^(\d+)-ISBL";
        press_sort_direction = "descend";
        metadata_columns = ["Element", "Unit", "Comment"];
        coordinate_order = ["t", "y", "x", "p", "q"];
    end
    % indexing - go into data() and cut time, raster, metadata, pressures accordingly

    properties (Dependent)
        grid_vectors (1, 5) cell;
        forecast (1, :) duration;
    end
    
%% PUBLIC OBJECT METHODS
    methods
        function obj = nwpdata(path, params)
            arguments
                path (1,:) string;
                params.elements (1,:) string;
                params.lat (1,1) double = 0;
                params.lon (1,1) double = 0;
                params.side (1,1) double = Inf;
                params.levels (1,:) = [-Inf Inf];
            end

            if path == ""
                obj.data = NaN(0, 0, 0, 0, 0);
                obj.time = datetime.empty;
                obj.raster = map.rasterref.MapPostingsReference.empty;
                obj.levels = [];
                obj.metadata = table;
                obj.created = datetime.empty;
                return; % get empty values
            end

            % Initialization
            elements = params.elements;
            levels = params.levels;
            
            info = georasterinfo(path);
            metadata = info.Metadata;
            raster = info.RasterReference;
            [obj.raster, x_idx, y_idx] = nwpdata.cropraster(raster, ...
                params.lat, params.lon, params.side);

            % YYYYY-ISBL data sampling
            if isnumeric(levels) 
                [all_pressures, ~] = nwpdata.find_press(metadata, elements, levels);
                obj.levels = unique(all_pressures, "stable");
                obj.data = NaN([1, obj.raster.RasterSize, ...
                    length(obj.levels), length(elements)]);
                % Now enough is known to allocate the object data member

                % read elements in one at a time - all-at-once takes too much memory
                for i_element = 1:length(elements)
                    [pressures, bands] = nwpdata.find_press(metadata, elements(i_element), levels);
                    [data, ~] = readgeoraster(path, Bands = bands);

                    % "Shouldn't" be necessary to intersect and match indices
                    % since the data appears to be consistently sorted and
                    % sampled at all levels, but just to be safe
                    [~, obj_idx, data_idx] = intersect(obj.levels, pressures, "stable");
                    if length(obj_idx) < size(obj.data, 4)
                        warning("Missing data for element %s at %d pressures", ...
                            elements(i_element), size(obj.data, 4) - length(obj_idx));
                    end

                    % Put data into object
                    obj.data(1, :, :, obj_idx, i_element) = data(y_idx, x_idx, data_idx);
                    obj.metadata = [obj.metadata; metadata(bands(1), :)];
                end
            else % Custom levels
                all_bands = nwpdata.find_bands(metadata, elements, levels);
                obj.levels = unique(metadata.ShortName(all_bands));
                obj.data = NaN([1, obj.raster.RasterSize, 1, length(all_bands)]);

                for i_element = 1:length(elements)
                    bands = nwpdata.find_bands(metadata, elements(i_element), levels);
                    this_levels = metadata.ShortName(bands);
                    [data, ~] = readgeoraster(path, Bands = bands);
                    [~, obj_idx, data_idx] = intersect(obj.levels, this_levels, "stable");
                    if length(obj_idx) < size(obj.data, 4)
                        warning("Missing data for element %s at %d levels", ...
                            elements(i_element), size(obj.data, 4) - length(obj_idx));
                    end
                    % Put data into object
                    obj.data(1, :, :, obj_idx, i_element) = data(y_idx, x_idx, data_idx);
                    obj.metadata = [obj.metadata; metadata(bands(1), :)];
                end
                % [data, ~] = readgeoraster(path, Bands = bands);
                %
                % obj.data(1, :, :, 1, :) = data(y_idx, x_idx, :);
                % obj.metadata = metadata(bands, :);
                % obj.levels = [];
            end
            obj.time = metadata.ValidTime(1);
            obj.created = metadata.ReferenceTime(1);
            obj.metadata = obj.metadata(:, ["Element", "Unit", "Comment"]);
        end

        function grid_vectors = get.grid_vectors(obj)
            [x, y] = worldGrid(obj.raster, "gridvectors");
            grid_vectors = {obj.time, y, x, obj.levels, obj.metadata.Element};
        end

        function forecast = get.forecast(obj)
            forecast = obj.time - obj.created;
        end

        function mapshow(obj)
            if any(size(obj.data, [1 4 5]) ~= 1)
                error("Data must be singleton across all non-xy dimensions");
            end
            map_data = squeeze(obj.data(1, :, :, 1, 1));
            
            mapshow(map_data, obj.raster, DisplayType = "surface");
            cb = colorbar;
            cb.Label.String = obj.metadata.Comment;
            axis equal;
            axis tight;
        end
    end

%% INDEXING METHODS
methods (Access=protected)
    % The acutal indexing logic lives in (index_vectors) and are further broken
    % out into index_(...) for maintainability
    function obj = parenReference(obj, indexOp)
        indices = indexOp.Indices;
        times = obj.time_index(indices{1});
        lvls = obj.lvl_index(indices{4});
        qties = obj.qty_index(indices{5});

        obj.data = obj.data(times, indices{2:3}, lvls, qties);
        obj.raster = obj.yx_index(indices{2:3});
        obj.time = obj.time(times);
        obj.created = obj.created(times);
        obj.levels = obj.levels(lvls);
        obj.metadata = obj.metadata(qties, :);
    end

    function obj = parenAssign(obj,indexOp,varargin)
        error("Paren assignment not supported");
        % % Ensure object instance is the first argument of call.
        % if isempty(obj)
        %     obj = varargin{1};
        % end
        % if isscalar(indexOp)
        %     assert(nargin==3);
        %     rhs = varargin{1};
        %     obj.ContainedArray.(indexOp) = rhs.ContainedArray;
        %     return;
        % end
        
        % [obj.(indexOp(2:end))] = varargin{:};
    end

    function n = parenListLength(~, ~, ~)
        n = 1;
    end

    function obj = parenDelete(~, ~, ~)
        error("Paren deletion not supported");
    end

    function time_indices = time_index(obj, time_op)
        switch (class(time_op))
            case "duration"
                [~, time_indices, ~] = intersect(time_op, obj.forecast, "stable");
            case "datetime"
                [~, time_indices, ~] = intersect(time_op, obj.time, "stable");
            case "timerange"
                tr_struct = struct(time_op);

                if isduration(tr_struct.first)
                    dummy = timetable((1:length(obj.time))', RowTimes = obj.forecast);
                elseif isdatetime(tr_struct.first)
                    dummy = timetable((1:length(obj.time))', RowTimes = obj.time);
                else
                    error("Unrecognized type for timerange subscript: %s", ...
                        class(tr_struct.first));
                end

                time_indices = dummy{time_op, 1};
            otherwise
                time_indices = time_op;
        end
    end

    function raster = yx_index(obj, y_op, x_op)
        grids = obj.grid_vectors;
        x = grids{3}(x_op);
        y = grids{2}(y_op);
        
        if isempty(x) || isempty(y)
            raster = map.rasterref.MapPostingsReference.empty;
            return;
        end
        
        x_limits = [min(x) max(x)];
        y_limits = [min(y) max(y)];
        
        raster = nwpdata.cropraster(obj.raster, x_limits, y_limits);
    end

    function lvl_indices = lvl_index(obj, lvl_op)
        if isstring(lvl_op)
            if isnumeric(obj.levels)
                error("String-based indexing not supported for numeric pressure levels");
            end
            [~, ~, lvl_indices] = intersect(lvl_op, obj.levels);
        else
            lvl_indices = lvl_op;
        end
    end

    function qty_indices = qty_index(obj, qty_op)
        if isstring(qty_op)
            [~, ~, qty_indices] = intersect(qty_op, obj.metadata.Element);
        else
            qty_indices = qty_op;
        end
    end
end 

methods (Access=public)
    function out = value(obj)
        out = squeeze(obj.data);
    end

    % function out = sum(obj)
    %     error("Sum not supported")
    %     % out = sum(obj.ContainedArray,"all");
    % end

    function out = cat(dim,varargin)
        if dim ~= 1
            error("NWP datasets can only be concatenated along the time dimension");
        end

        out = timecat(varargin{:});
    end

    function out = horzcat(varargin)
        out = timecat(varargin{:});
    end

    function obj = timecat(varargin)
        for i = 2:length(varargin)
            ref = varargin{1};
            objut = varargin{i};
            if ~nwpdata.rasterequals(ref.raster, objut.raster)
                error("Raster mismatch at position %d", i);
            end
            if any(ref.levels ~= objut.levels)
                error("Level mismatch at position %d");
            end
            if any(ref.metadata.Element ~= objut.metadata.Element)
                error("Quantity mismatch at position %d");
            end
        end

        datas = cellfun(@(obj) obj.data, varargin, UniformOutput = false);
        times = cellfun(@(obj) obj.time, varargin, UniformOutput = false);
        obj = varargin{1};
        obj.data = cat(1, datas{:});
        obj.time = cat(2, times{:});
    end

    function varargout = size(obj,varargin)
        [varargout{1:nargout}] = size(obj.data,varargin{:});
    end
end

methods (Static, Access=public)
    function obj = empty()
        obj = nwpdata("");
    end
end
    
%% PUBLIC UTILITIES
methods (Static)
        function [raster, x_idx, y_idx] = cropraster(raster, varargin)
            % [raster, x_idx, y_idx] = CROPRASTER(raster, lat, lon, side)
            % Crop the (raster) object to a square centered at (lat, lon) that is (side) wide
            % Return the logical index vectors to crop the associated matrix
            % (separated out because they are used multiple times)

            [x_values, y_values] = worldGrid(raster, "gridvectors");
            switch nargin
                case 3
                    x_limits = varargin{1};
                    y_limits = varargin{2};
                case 4
                    lat = varargin{1};
                    lon = varargin{2};
                    side = varargin{3};
                    if side == Inf
                        x_idx = true(1, raster.RasterSize(2));
                        y_idx = true(1, raster.RasterSize(1));
                        return;
                    end

                    [x_center, y_center] = projfwd(raster.ProjectedCRS, lat, lon);
                    x_limits = x_center + side/2*[-1 1];
                    y_limits = y_center + side/2*[-1 1];
                otherwise
                    error("Unrecognzied number of input arguments");
            end

            [~, raster] = mapcrop(zeros(raster.RasterSize), raster, x_limits, y_limits);
            if isempty(raster)
                error("Cropping produced empty result");
            end

            x_idx = raster.XWorldLimits(1) <= x_values & x_values <= raster.XWorldLimits(2);
            y_idx = raster.YWorldLimits(1) <= y_values & y_values <= raster.YWorldLimits(2);

            assert(sum(x_idx) == raster.RasterSize(2), "X indexor must match raster size");
            assert(sum(y_idx) == raster.RasterSize(1), "Y indexor must match raster size");
        end

        % function [raster, x_offset, y_offset] = centerraster(raster, lat, lon)
        %     arguments
        %         raster
        %         lat (1,1) double = NaN;
        %         lon (1,1) double = NaN;
        %     end 
        % end

        function isequal = rasterequals(raster1, raster2)
            match_xlimits = all(raster1.XWorldLimits == raster2.XWorldLimits);
            match_ylimits = all(raster1.YWorldLimits == raster2.YWorldLimits);
            match_size = all(raster1.RasterSize == raster2.RasterSize);
            match_crs = raster1.ProjectedCRS.isequal(raster2.ProjectedCRS);
            isequal = match_xlimits && match_ylimits && match_size && match_crs;
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

        function [pressures, bands] = find_press(metadata, elements, levels)
            % [pressures, bands] = FIND_PRESS(metadata, elements, levels)
            % Find the unique pressure levels associated with (elements) isolated to range (levels)
            % Sort in the order configured under Constants
            arguments
                metadata table;
                elements (1,:) string;
                levels (1,2) double;
            end
            bands = nwpdata.find_bands(metadata, elements, nwpdata.press_regex);
            pressures = string(regexp(metadata.ShortName(bands), nwpdata.press_regex, "tokens"));
            pressures = str2double(pressures);
            % cut levels to specified range
            in_range = min(levels) <= pressures & pressures <= max(levels);
            pressures = pressures(in_range);
            bands = bands(in_range);

            [pressures, order] = sort(pressures, nwpdata.press_sort_direction);
            bands = bands(order);
        end

        % Get the URL for a NAM analysis
        function [filename, subfolder, url] = make_url(date)
            
            base = "https://www.ncei.noaa.gov/data/north-american-mesoscale-model/access/analysis/";
            folder_fmt = "yyyyMM/yyyyMMdd/";
            file_fmt = "'nam_218'_yyyyMMdd_HHmm_000.'grb2'";
            subfolder = string(date, folder_fmt);
            filename = string(date, file_fmt);
            url = base + subfolder + filename;
        end
        %
        % function [files, downloaded] = download(dates, dest)
        %     [names, folders, urls] = make_url(dates);
        %     for i = 1:length(urls)
        %         
        %     end
        % end
    end

%% INTERNAL UTILITIES
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
