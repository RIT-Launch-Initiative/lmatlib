% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

% Methods of indexing:
% 

classdef (Sealed) nwpdata < matlab.mixin.indexing.RedefinesParen

    properties (Dependent, Access = public)
        data double;
        time (1, :) duration;
        proj_x (1, :) double;
        proj_y (1, :) double;
        layers (1, :) string;
        quantities (1, :) string;
        epoch (1, 1) datetime;
        projection (1, 1) projcrs;
    end

    properties (Access = private)
        DATA (:, :, :, :, :) double; 
        PROJ_X (1, :) double;
        PROJ_Y (1, :) double;
        DATETIME (1, :) datetime;
        LEVELS (1, :) string;
        METADATA table;
        PROJECTION (1, 1) projcrs;
    end

    properties (Constant, Access = private)
        metadata_columns = ["Description", "Element", "Unit", "Comment"];
        coordinate_order = ["time", "proj_x", "proj_y", "layer", "qty"];
        grib_coordinate_order = ["proj_y", "proj_x", "layer"];
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
                params.layers (1,:) string;
            end

            if path == ""
                return; % get empty values
            end

            % Initialization
            elements = params.elements;
            layers = params.layers;
            
            info = georasterinfo(path);
            metadata = info.Metadata;
            raster = info.RasterReference;
            obj.PROJECTION = raster.ProjectedCRS;
            [raster, x_idx, y_idx] = nwpdata.cropraster(raster, ...
                params.lat, params.lon, params.side);
            [obj.PROJ_X, obj.PROJ_Y] = worldGrid(raster, "gridvectors");

            all_bands = nwpdata.find_bands(metadata, elements, layers);
            obj.LEVELS = unique(metadata.ShortName(all_bands));

            % NOTE assumes indices (t, x, y, l, q)
            obj.DATA = NaN([1, length(obj.PROJ_X), length(obj.PROJ_Y), ...
                length(obj.LEVELS), length(elements)]);

            for i_element = 1:length(elements)
                bands = nwpdata.find_bands(metadata, elements(i_element), layers);
                this_layers = metadata.ShortName(bands);
                [data, ~] = readgeoraster(path, Bands = bands);

                [~, obj_idx, data_idx] = intersect(obj.LEVELS, this_layers, "stable");
                if length(obj_idx) < length(obj.LEVELS)
                    warning("Missing %s data at %s", ...
                        elements(i_element), setdiff(obj.LEVELS, this_layers, "stable"));
                end
                
                [~, ~, data_reorder] = intersect(nwpdata.coordinate_order, ...
                    nwpdata.grib_coordinate_order, "stable");

                data = data(y_idx, x_idx, data_idx);
                data = permute(data, data_reorder);

                % Put data into object
                
                indices = nwpdata.ordered_index(time = 1, layer = obj_idx, qty = i_element);
                obj.DATA(indices{:}) = data;
                obj.METADATA = [obj.METADATA; metadata(bands(1), nwpdata.metadata_columns)];
            end

            obj.DATETIME = metadata.ValidTime(1);
        end

        function mapshow(obj)
            [~, ~, obj_reorder] = intersect(nwpdata.grib_coordinate_order, ...
                nwpdata.coordinate_order, "stable");

            map_data = permute(squeeze(obj.DATA), obj_reorder);
            if ~ismatrix(map_data)
                error("Too many non_singleton dimensions (%d) to plot (%d)", ndims(map_data), 2);
            end

            mapshow(obj.proj_x, obj.proj_y, map_data, DisplayType = "surface");
            cb = colorbar;
            cb.Label.String = obj.METADATA.Comment;
            axis equal;
            axis tight;
        end
    end

%% GETTERS
    methods
        function data = get.data(obj)
            data = squeeze(obj.DATA);
        end

        function x = get.proj_x(obj)
            x = obj.PROJ_X;
        end

        function y = get.proj_y(obj)
            y = obj.PROJ_Y;
        end

        function t = get.epoch(obj)
            t = obj.DATETIME(1);
        end

        function t = get.time(obj)
            t = obj.DATETIME - obj.epoch;
        end

        function lev = get.layers(obj)
            lev = obj.LEVELS;
        end

        function crs = get.projection(obj)
            crs = obj.PROJECTION;
        end

        function quants = get.quantities(obj)
            quants = obj.METADATA.Element';
        end
    end

%% INDEXING METHODS
methods (Access=protected)

    function obj = parenReference(obj, indexOp)
        indices = indexOp.Indices;

        % indices = struct("indices");

        time_indices = obj.time_index(indices{1});
        x_indices = indices{2};
        y_indices = indices{3};
        layer_indices = obj.layer_index(indices{4});
        qty_indices = obj.qty_index(indices{5});
        data_indices = nwpdata.ordered_index(time = time_indices, ...
            proj_x = x_indices, proj_y = y_indices, ...
            layer = layer_indices, qty = qty_indices);

        obj.DATETIME = obj.DATETIME(time_indices);
        obj.PROJ_X = obj.PROJ_X(x_indices);
        obj.PROJ_Y = obj.PROJ_Y(y_indices);
        obj.LEVELS = obj.LEVELS(layer_indices);
        obj.METADATA = obj.METADATA(qty_indices, :);
        obj.DATA = obj.DATA(data_indices{:});
    end

    function parenAssign(~, ~, ~)
        error("Paren assignment not supported");
    end

    function n = parenListLength(~, ~, ~)
        n = 1;
    end

    function parenDelete(~, ~, ~)
        error("Paren deletion not supported");
    end

    function time_indices = time_index(obj, time_op)
        dummy = timetable((1:length(obj.time))', RowTimes = obj.time);
        time_indices = dummy{time_op, 1};
    end

    function layer_indices = layer_index(obj, layer_op)
        if isstring(layer_op)
            [~, ~, layer_indices] = intersect(layer_op, obj.LEVELS, "stable");
        else
            layer_indices = layer_op;
        end
    end

    function qty_indices = qty_index(obj, qty_op)
        if isstring(qty_op)
            [~, ~, qty_indices] = intersect(qty_op, obj.METADATA.Element, "stable");
        else
            qty_indices = qty_op;
        end
    end

end 

methods (Access=public)
    function out = value(obj)
        out = obj.DATA;
    end

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
        nwpdata.compatible(varargin{:});

        datas = cellfun(@(obj) obj.DATA, varargin, UniformOutput = false);
        times = cellfun(@(obj) obj.DATETIME, varargin, UniformOutput = false);
        obj = varargin{1};
        obj.DATA = cat(find(nwpdata.coordinate_order == "time"), datas{:});
        obj.DATETIME = [times{:}];

        if ~issorted(obj.DATETIME, "strictascend")
            error("Time must strictly increase");
        end
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

        function bands = find_bands(metadata, elements, layers)
            % [bands] = FIND_BANDS(metadata, elements, layers)
            % Find the bands (linear indices) matching any elements AND any layers
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
            bands = find(element_idx & layer_idx);
        end

        function ops = ordered_index(ops)
            arguments
                ops.time = ':';
                ops.proj_x = ':';
                ops.proj_y = ':';
                ops.layer = ':';
                ops.qty = ':';
            end
            
            ops = arrayfun(@(field) ops.(field), nwpdata.coordinate_order, ...
                UniformOutput = false);
        end

        % Get the URL for a NAM analysis
        function [filename, subfolder, url] = make_url(date)
            if any(mod(date.Hour, 6))
                error("Invalid model hour");
            end

            base = "https://www.ncei.noaa.gov/data/north-american-mesoscale-model/access/analysis/";
            folder_fmt = "yyyyMM/yyyyMMdd/";
            file_fmt = "'nam_218'_yyyyMMdd_HHmm_000.'grb2'";
            subfolder = string(date, folder_fmt);
            filename = string(date, file_fmt);
            url = base + subfolder + filename;
        end
        
        function [files, downloaded] = download(dates, dest)
            [names, ~, urls] = nwpdata.make_url(dates);
            files = fullfile(dest, names);
            downloaded = false(size(files));

            for i = 1:length(urls)
                fprintf("Downloading (%d) %s from %s\n", i, names(i), urls(i));
                if isfile(files(i))
                    downloaded(i) = true;
                    fprintf("Located on disk at %s\n", files(i));
                else
                    tic;
                    downloaded(i) = copyfile(urls(i), files(i));
                    timed = toc;
                    fprintf("Finished in %.2f sec\n", timed);
                end

                if ~downloaded(i)
                    warning("Unable to download '%s' at '%s'", names(i), urls(i));
                end
            end
            files = files(downloaded);
        end
    end

%% INTERNAL UTILITIES
    methods (Static, Access = public)
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

        function compatible(varargin)
            err_id = "nwpdata:cat";
            for i = 2:length(varargin)
                ref = varargin{1};
                objut = varargin{i};
                tests = [isequal(ref.proj_x, objut.proj_x), ...
                    isequal(ref.proj_y, objut.proj_y), ...
                    isequal(ref.layers, objut.layers), ...
                    isequal(ref.quantities, objut.quantities)];
                if any(~tests)
                    mex = MException(err_id, "Dimension mismatch at %d", i);
                    throwAsCaller(mex)
                end
            end
        end
    end

end
