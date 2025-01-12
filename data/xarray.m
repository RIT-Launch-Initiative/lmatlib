% Array with named coordinate variables
% MATLAB implementation of xarray
% Paren indexing: index into data and underlying axisinfo "normally"
% Brace indexing: index into data only
% Dot indexing: index into axes

classdef xarray < matlab.mixin.indexing.RedefinesDot ...
        & matlab.mixin.indexing.RedefinesParen ...
        & matlab.mixin.indexing.RedefinesBrace ...
        & matlab.mixin.CustomDisplay

    properties (GetAccess = public, SetAccess = protected)
        axes (1, :) string = [];
        coordinates (1, :) cell = {};
    end

    properties (Access = protected)
        data double = [];
    end

    methods
        function obj = xarray(data, axes, coordinates)
            arguments
                data double;
            end

            arguments (Repeating)
                axes (1, :) string;
                coordinates (1, :);
            end
            
            if ~isscalar(axes{1}) && iscell(coordinates{1})
                axes = axes{1};
                coordinates = coordinates{1};
            end

            if length(axes) < ndims(data)
                error("Not enough coordinate axes specified");
            end

            obj.data = data;
            for i_ax = 1:length(axes)
                obj.(axes{i_ax}) = coordinates{i_ax};
            end
        end

        function data = double(obj)
            data = squeeze(obj.data);
        end
    end

%% MANIUPLATION
    methods (Access = public)
        function obj = cat(dim, arrays)
            arguments
                dim (1,1) {mustBeA(dim, ["numeric", "string"])};
            end
            arguments (Repeating)
                arrays xarray;
            end
            template = arrays{1};
            dim = template.convertdims(dim);

            to_compare = 1:naxes(template);
            to_compare(dim) = [];

            for i_obj = 2:length(arrays)
                if ~dimequal(template, arrays{i_obj}, to_compare)
                    error("Mismatched dimensions at position %d", i_obj);
                end
            end

            datas = cellfun(@(arry) arry.data, arrays, UniformOutput = false);
            coords = cellfun(@(arry) arry.coordinates{dim}, arrays, UniformOutput = false);
            template.data = cat(dim, datas{:});
            template.coordinates{dim} = vertcat(coords{:});
            obj = template;
        end


        function obj = permute(obj, order)
            dimorder = obj.convertdims(order);

            if length(dimorder) ~= naxes(obj)
                error("Dimension order must have exactly one entry per axis:\n%s", ...
                    mat2str(obj.axes));
            end

            obj.data = permute(obj.data, dimorder);
            obj.axes = obj.axes(dimorder);
            obj.coordinates = obj.coordinates(dimorder);
        end

        function obj = squeeze(obj)
            scalardims = size(obj) == 1;
            obj.data = squeeze(obj.data);
            obj.axes(scalardims) = [];
            obj.coordinates(scalardims) = [];
        end

        % Index by exact equality, according to ismember(<values>, <coord>)
        % Returns values in the order specified by <values>
        % xarr.pick(axis = "value")
        function obj = pick(obj, axes, values)
            arguments
                obj xarray
            end
            arguments (Repeating)
                axes (1, 1) string;
                values (:, 1);
            end

            dims = obj.convertdims(axes{:});
            ops = repmat({':'}, 1, naxes(obj));

            for i_ax = 1:length(axes)
                [ops{dims(i_ax)}, ~] = find(values{i_ax}' == obj.coordinates{i_ax});
                % [~, ops{dims(i_ax)}] = ismember(values{i_ax}, obj.(axes{i_ax}));
                if ops{dims(i_ax)} == 0
                    ops{dims(i_ax)} = [];
                end
            end

            obj = obj(ops{:});
        end

        % Index by equality with tolerance
        % xarr.range(axis = [center tol])
        function obj = pickt(obj, axes, ranges)
            arguments
                obj xarray
            end
            arguments (Repeating)
                axes (1, 1) string;
                ranges (1, 2);
            end
            
            ops = repmat({':'}, 1, naxes(obj));
            dims = obj.convertdims(axes{:});

            for i_ax = 1:length(axes)
                coord = obj.(axes{i_ax});
                rang = ranges{i_ax}(1) + ranges{i_ax}(2) * [-1 1];
                ops{dims(i_ax)} = rang(1) <= coord & coord <= rang(2);
            end

            obj = obj(ops{:});
        end

        % Index by range membership, according to min(rng) <= coord & coord <= max(rng)
        % xarr.range(axis = [lower, upper])
        %   [upper, lower] also works, and either can be +/- Inf
        function obj = range(obj, axes, ranges)
            arguments
                obj xarray
            end
            arguments (Repeating)
                axes (1, 1) string;
                ranges (1, 2);
            end

            ops = repmat({':'}, 1, naxes(obj));
            dims = obj.convertdims(axes{:});

            for i_ax = 1:length(axes)
                coord = obj.(axes{i_ax});
                ops{dims(i_ax)} = min(ranges{i_ax}) <= coord & coord <= max(ranges{i_ax});
            end

            obj = obj(ops{:});
        end


        function obj = index(obj, axes, indices)
            arguments
                obj xarray
            end
            arguments (Repeating)
                axes (1, 1) string;
                indices (1, :);
            end

            ops = repmat({':'}, 1, naxes(obj));
            dims = obj.convertdims(axes{:});

            for i_ax = 1:length(axes)
                ops{dims(i_ax)} = indices{i_ax};
            end

            obj = obj(ops{:});
        end
    end


%% METADATA ACCESS
    methods (Access = public)
        function varargout = size(obj, varargin)
            if isempty(varargin)
                sizes = cellfun(@length, obj.coordinates);
            else
                queries = obj.convertdims(varargin{:});
                sizes = cellfun(@length, obj.coordinates(queries));
            end
            
            varargout = cell(1, nargout);
            if nargout == 1
                varargout{1} = sizes;
            else
                for i = 1:nargout
                    varargout{i} = sizes(i);
                end
            end
        end

        function n = naxes(obj)
            n = length(obj.axes);
        end

        function obj = sort(obj, dims, varargin)
            dims = obj.convertdims(dims);

            ops = repmat({':'}, 1, naxes(obj));
            for dim = dims
                [~, ops{dim}] = sort(obj.coordinates{dim}, varargin{:});
            end
            obj = obj(ops{:});
        end
    end

    methods (Static)
        function obj = empty
            obj = xarray([]);
        end
    end

    methods (Access = protected)
%% INDEXING HELPER FUNCTIONS
        function dims = convertdims(obj, varargin)
            % array from input
            dims = [varargin{:}];
            if isnumeric(dims)
                return;
            end

            % dims() is a string
            names = dims;
            matches = dims == obj.axes';
            [dims, ~] = find(matches);
            if length(dims) < length(names)
                mex = MException("xarray:invalidAxis", ...
                    "Unrecognized axis or property name(s) '%s'", names(~any(matches, 1)));
                throwAsCaller(mex);
            end
        end

        function [iseq] = dimequal(obj, other, dims)
            iseq = false(1, length(dims));
            for dim = dims
                iseq(dim) = obj.axes(dim) == other.axes(dim) &&...
                    isequal(obj.coordinates{dim}, other.coordinates{dim});
            end
        end

        function indexcheck(obj, op)
            if length(op.Indices) ~= naxes(obj)
                mex = MException("xarray:indexing", ...
                    "%d subscripts expected, got %d", ...
                    naxes(obj), length(op.Indices));
                throwAsCaller(mex);
            end
        end

%% INDEXING OVERRIDES
        function obj = parenReference(obj, operations)
            % Input processing and checking
            op = operations(1);
            indexcheck(obj, op);

            % Index into data using all indices, index into each axis using its index
            obj.data = obj.data.(op);
            for i_axis = 1:naxes(obj)
                obj.coordinates{i_axis} = obj.coordinates{i_axis}(op.Indices{i_axis});
            end

            % Perform subsequent operations (if any)
            if ~isscalar(operations)
                obj = obj.(operations(2:end));
            end
        end

        function obj = parenAssign(obj, operations, subobj)
            % Input processing and checking
            if ~isscalar(operations)
                error("Chained paren-assignment not supported")
            end
            op = operations(1);
            indexcheck(obj, op);

            % Test sub-object coordinates for equality
            for i_axis = 1:naxes(obj)
                coord = obj.coordinates{i_axis};
                if ~isequal(coord(op.Indices{i_axis}), subobj.coordinates{i_axis})
                    error("Unable to assign sub-array: mismatch along %s", obj.names(i_axis));
                end
            end

            % Apply operation
            obj.data.(op) = subobj.data;
        end

        function obj = parenDelete(obj, operations)
            % Input processing and checking
            if ~isscalar(operations)
                error("Chained paren-deletion not supported")
            end
            op = operations(1);
            indexcheck(obj, op);

            for i_axis = 1:naxes(obj)
                obj.coordinates{i_axis}(op.Indices{i_axis}) = [];
            end

            obj.data.(op) = [];
        end
        
        function n = parenListLength(~, ~, ~)
            n = 1;
        end

        function data = braceReference(obj, operations)
            % Input processing and checking
            if ~isscalar(operations)
                error("Chained brace-reference not supported")
            end

            op = operations(1);
            data = obj.data(op.Indices{:});
        end
        
        function obj = braceAssign(obj, op, data)
            obj.data(op.Indices{:}) = data;
        end

        function n = braceListLength(~, ~, ~)
            n = 1;
        end

        % OVERRIDE DOT-INDEXING TO ACCESS ARRAY AXES
        function coord = dotReference(obj, operations)
            op = operations(1);

            dim = obj.convertdims(op.Name);
            coord = obj.coordinates{dim};

            if ~isscalar(operations)
                coord = coord.(operations(2:end));
            end
        end

        function obj = dotAssign(obj, op, input)
            arguments
                obj xarray;
                op (1,1) matlab.indexing.IndexingOperation;
                input (:, 1);
            end

            axis = op.Name;


            axis_idx = find(axis == obj.axes);
            if isempty(axis_idx)
                if any(axis == properties(obj)) || any(axis == methods(obj))
                    error("'%s' already used as object property or method", axis);
                end
                axis_idx = length(obj.axes) + 1; % create a new axis if it doesn't exist
            end

            if length(input) ~= size(obj.data, axis_idx)
                error("Dimension mismatch: size(data, %d) = %d, length(%s) = %d", ...
                    axis_idx, size(obj.data, axis_idx), axis, length(input));
            end

            obj.axes(axis_idx) = axis;
            obj.coordinates{axis_idx} = input;
        end

        function n = dotListLength(~, ~, ~)
            n = 1;
        end
    end

%% OBJECT DISPLAY OVERRIDE
    methods (Access = protected)
        function header = getHeader(obj)
            dimstr = matlab.mixin.CustomDisplay.convertDimensionsToString(obj);
            nonscalar = sum(size(obj) ~= 1);
            if naxes(obj) > 3 && nonscalar ~= naxes(obj)
                dimstr = sprintf('%s (%d nonscalar)', dimstr, nonscalar);
            end
            header = sprintf('%s %s with axes:', dimstr, class(obj));
        end

        function group = getPropertyGroups(obj)
            prop_list = [cellstr(obj.axes); 
                cellfun(@(c) c', obj.coordinates, UniformOutput = false)];
            prop_list = struct(prop_list{:});
            group = matlab.mixin.util.PropertyGroup(prop_list);
        end
    end
end

