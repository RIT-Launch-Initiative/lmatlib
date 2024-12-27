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
                data {mustBeA(data, ["double", "xarray"])};
            end

            arguments (Repeating)
                axes (1, 1) string;
                coordinates (1, :);
            end
            
            if isa(data, "xarray")
                obj = data;
                return;
            end

            if isempty(data)
                return;
            end
            if length(axes) < ndims(data)
                error("Not enough coordinate axes specified");
            end

            obj.data = data;
            for i_ax = 1:length(axes)
                obj.(axes{i_ax}) = coordinates{i_ax};
            end
        end

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

        function obj = cat(dim, varargin)
            obj = varargin{1};
            dim = convertdims(obj, dim);

            to_compare = 1:ndims(varargin{1});
            to_compare = to_compare(to_compare ~= dim);
            for i_obj = 2:length(varargin)
                if ~obj.dimequal(varargin{i_obj}, to_compare)
                    error("Mismatched dimensions at position %d", i_obj);
                end
            end

            datas = cellfun(@(obj) obj.data, varargin, UniformOutput = false);
            coords = cellfun(@(obj) obj.coordinates{dim}, varargin, UniformOutput = false);
            obj.data = cat(dim, datas{:});
            obj.coordinates{dim} = [coords{:}];
        end

        function n = ndims(obj)
            n = length(obj.axes);
        end

        function obj = sort(obj, dims, varargin)
            dims = obj.convertdims(dims);

            ops = repmat({':'}, 1, ndims(obj));
            for dim = dims
                [~, ops{dim}] = sort(obj.coordinates{dim}, varargin{:});
            end
            obj = obj(ops{:});
        end

        function obj = permute(obj, dimorder)
            dimorder = obj.convertdims(dimorder);

            if length(dimorder) ~= ndims(obj)
                error("Dimension order must have exactly one entry per coordinate");
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

        function data = double(obj)
            data = squeeze(obj.data);
        end

        function obj = subset(obj, axes, indices)
            arguments
                obj xarray
            end
            arguments (Repeating)
                axes (1, 1) string;
                indices (1, :);
            end

            ops = repmat({':'}, 1, ndims(obj));
            dims = obj.convertdims(axes{:});

            for i_ax = 1:length(axes)
                if isnumeric(indices{i_ax}) || islogical(indices{i_ax})
                    ops{dims(i_ax)} = indices{i_ax};
                else
                    [~, ops{dims(i_ax)}] = ismember(indices{i_ax}, obj.(axes{i_ax}));
                end
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
        % OVERRIDE PAREN REFERENCING TO SUBSCRIPT AXES
        function dims = convertdims(obj, varargin)
            dims = [varargin{:}];
            if isstring(dims)
                names = dims;
                [present, dims] = ismember(names, obj.axes);
                if any(~present)
                    error("Unrecognized axes: %s", mat2str(names(~present)));
                end
            end
        end

        function [iseq] = dimequal(obj, other, dims)
            iseq = false(1, length(dims));
            for dim = dims
                iseq(dim) = obj.axes(dim) == other.axes(dim) &&...
                    isequal(obj.coordinates{dim}, other.coordinates{dim});
            end
        end

        function obj = parenReference(obj, op)
            if length(op.Indices) ~= ndims(obj)
                error("%d subscripts expected, got %d", ndims(obj), length(op.Indices));
            end

            obj.data = obj.data.(op);
            for i = 1:ndims(obj)
                obj.(obj.axes(i)) = obj.(obj.axes(i))(op.Indices{i});
            end
        end

        function obj = parenAssign(obj, op, subobj)
            if length(op.Indices) ~= ndims(obj)
                error("%d subscripts expected, got %d", ndims(obj), length(op.Indices));
            end

            dest = obj.(op);
            if ~all(dest.dimequal(subobj, 1:ndims(dest)))
                error("Dimensions not equal");
            end
            obj.data(op.Indices{:}) = subobj.data;
        end

        function parenDelete(~, ~)
            error("xarray does not support paren-deletion");
        end
        
        function n = parenListLength(~, ~, ~)
            n = 1;
        end

        function data = braceReference(obj, op)
            data = obj.data(op.Indices{:});
        end
        
        function obj = braceAssign(obj, op, data)
            obj.data(op.Indices{:}) = data;
        end

        function n = braceListLength(~, ~, ~)
            n = 1;
        end

        % OVERRIDE DOT-INDEXING TO ACCESS ARRAY AXES
        function coord = dotReference(obj, op)
            axis = op(1).Name;
            axis_idx = find(axis == obj.axes);
            if isempty(axis_idx)
                error("No axis '%s'", axis);
            end
            % if the data is subsequently time-indexed
            if length(op) == 2
                coord = obj.coordinates{axis_idx}.(op(2));
            else
                coord = obj.coordinates{axis_idx};
            end
        end

        function obj = dotAssign(obj, op, input)
            arguments
                obj xarray;
                op matlab.indexing.IndexingOperation;
                input (1, :);
            end
            axis = op.Name;
            if ismember(axis, ["data", "axes", "coordinates"])
                error("'%s' already used as xarray property name", axis);
            end
            if ~isvector(input)
                error("Coordinates must be a vector");
            end

            axis_idx = find(axis == obj.axes);
            if isempty(axis_idx)
                axis_idx = length(obj.axes) + 1;
            end

            if length(input) ~= size(obj.data, axis_idx)
                error("Mismatch in dimension %d: Data is %d, coordinates are %d", ...
                    axis_idx, size(obj.data, axis_idx), length(input));
            end

            obj.axes(axis_idx) = axis;
            obj.coordinates{axis_idx} = input;
        end

        function n = dotListLength(~, ~, ~)
            n = 1;
        end
    end

    methods (Access = protected)
        function header = getHeader(obj)
            dimstr = matlab.mixin.CustomDisplay.convertDimensionsToString(obj);
            nonscalar = sum(size(obj) ~= 1);
            if ndims(obj) > 3 && nonscalar ~= ndims(obj)
                dimstr = sprintf('%s (%d nonscalar)', dimstr, nonscalar);
            end
            header = sprintf('%s %s with axes:', dimstr, class(obj));
        end

        function group = getPropertyGroups(obj)
            prop_list = [cellstr(obj.axes); obj.coordinates];
            prop_list = struct(prop_list{:});
            group = matlab.mixin.util.PropertyGroup(prop_list);
        end
    end
end

