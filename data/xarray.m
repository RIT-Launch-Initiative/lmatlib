%% Array with named coordinate variables
% Barebones MATLAB implementation of python-xarray

classdef xarray < matlab.mixin.indexing.RedefinesDot ...
        & matlab.mixin.indexing.RedefinesParen ...
        & matlab.mixin.indexing.RedefinesBrace ...
        & matlab.mixin.CustomDisplay

    properties (GetAccess = public, SetAccess = protected)
        axes (1, :) string = [];
        coordinates (1, :) cell = {};
        data double = [];
    end

    methods
        function obj = xarray(data, axes, coordinates)
            % Construct array with named coordinate variables
            % xarr = xarray(data, axis1 = coordinates1, axis2 = coordinates2[, ...])
            %   data                double array
            %   axis1 ...           string axis name
            %   coordinates1 ...    coordinate variables (vectors, any type)
            % 
            % xarr = xarray(data, axes, coordinates)
            %   data                double array
            %   axes ...            string array of axis names
            %   coordinates ...     cell array of coordinate variables
            %
            % Coordinate variables must
            %   - Be vectors
            %   - Match the size of the data in that dimension
            arguments
                data double;
            end

            arguments (Repeating)
                axes (1, :) string;
                coordinates (1, :);
            end

            if nargin == 0
                obj.data = [];
                return;
            end
            
            if ~isscalar(axes{1}) && iscell(coordinates{1})
                axes = axes{1};
                coordinates = coordinates{1};
            end

            if length(axes) ~= length(coordinates)
                error("Specified %d axes but %d coordinates", ...
                    length(axes), length(coordinates));
            end

            if isvector(data) && (size(data, 2) > 1)
                data = data';
            end

            nd = ndims(data);
            if (nd == 2) && (size(data, 2) == 1)
                nd = 1;
            end

            if length(axes) < nd
                error("At least %d axes required but %d specified", ...
                    nd, length(axes))
            end

            obj.data = data;
            for i_axis = 1:length(axes)
                obj.(axes{i_axis}) = coordinates{i_axis};
            end
        end

        function dbl = double(obj)
            dbl = squeeze(obj.data);
        end

        function nax = naxes(obj)
            nax = length(obj.axes);
        end

        function varargout = size(obj, varargin)
            % Get the xarray's size by name or number
            % sz = size(xarr)
            %   get size along all dimensions
            % sz = size(xarr, [dim1, dim2, ...])
            % sz = size(xarr, dim1, dim2, ...)
            %   dim1...N can be strings or numbers, but must all be the same
            % 
            % Outputs can be produced in one vector or split into several
            % sz = size(xarr, __)
            % [sz1, sz2, ...] = size(xarr, __)
            if nargin == 1
                [varargout{1:nargout}] = size(obj.data);
            else
                dims = obj.convertdims(varargin{:});
                [varargout{1:nargout}] = size(obj.data, dims);
            end
        end
    end

%% MANIUPLATION
    methods (Access = public)
        function obj = cat(dim, arrays)
            % Concatenate xarray instances by named dimension
            % xarr = cat(dimension, array1, array2, ...)
            %   dimension   index or name
            arguments
                dim (1,1) {mustBeA(dim, ["numeric", "string"])};
            end
            arguments (Repeating)
                arrays xarray;
            end

            template = arrays{1};

            % process input arguments
            if isnumeric(dim)
                if dim > naxes(template)
                    error("Concatenation dimension %d too large (%d-D)." + ...
                        " To create a new dimension, specify the new axis as a string scalar")
                end
            else
                name = dim;
                dim = find(template.axes == name);

                if isempty(dim)
                    dim = naxes(template) + 1;
                end
            end
            
            to_compare = 1:naxes(template);

            % create new coordinate axis
            if dim <= naxes(template) % existing dimension
                new_axes = template.axes;
                coords = cellfun(@(arry) arry.coordinates{dim}, ...
                    arrays, UniformOutput = false);
                coords = vertcat(coords{:});
                to_compare(dim) = [];
            else % new dimension defaults to 1:n
                new_axes = [template.axes name];
                coords = 1:length(arrays);
            end


            for i_obj = 2:length(arrays)
                if any(~dimequal(template, arrays{i_obj}, to_compare))
                    error("Mismatched dimensions at position %d", i_obj);
                end
            end

            datas = cellfun(@(arry) arry.data, arrays, UniformOutput = false);
            template.data = cat(dim, datas{:});
            template.axes = new_axes;
            template.coordinates{dim} = coords;
            obj = template;
        end


        function obj = permute(obj, order)
            % Re-order axes by name or number
            % xarr = permute(xarr, [dim1 dim2 ...]);

            dimorder = obj.convertdims(order);

            if length(dimorder) ~= naxes(obj)
                error("Dimension order must have exactly one entry per axis:\n%s" + ...
                    "\nTo create singleton dimensions for broadcasting, assign 1-element axes.", ...
                    mat2str(obj.axes));
            end
            if issorted(dimorder)
                % Dimensions are already in order, do nothing
                return;
            end

            obj.data = permute(obj.data, dimorder);
            obj.axes = obj.axes(dimorder);
            obj.coordinates = obj.coordinates(dimorder);
        end

        function obj = leading(obj, order)
            % Turn specified dimensions into leading dimensions
            leading_order = obj.convertdims(order)';
            trailing_order = 1:naxes(obj);
            trailing_order(leading_order) = [];

            obj = permute(obj, [leading_order trailing_order]);
        end

        function obj = squeeze(obj)
            lengths = cellfun(@length, obj.coordinates);
            scalars = lengths == 1;
            obj.data = squeeze(obj.data);
            if isvector(obj.data) && (size(obj.data, 2) > 1)
                obj.data = obj.data';
            end
            obj.axes(scalars) = [];
            obj.coordinates(scalars) = [];
        end

        function obj = repmat(obj, axes, counts)
            arguments
                obj xarray;
            end
            arguments (Repeating)
                axes (1,1) string;
                counts (1,1) double {mustBeInteger};
            end

            axes = string(axes);
            repdims = obj.convertdims(axes);
            repcounts = ones(1, naxes(obj));
            repcounts(repdims) = [counts{:}];
            obj.data = repmat(obj.data, repcounts);
            for i_rep = 1:length(repdims)
                obj.coordinates{i_rep} = repmat(obj.coordinates{i_rep}, counts{i_rep}, 1);
            end
        end

        function [obj, indices] = pick(obj, axes, values)
            % Index by exact equality, according to ismember(<values>, <coord>)
            % Returns values in the order specified by <values>
            % xarr.pick(axis = "value", ...)
            
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
                [~, ops{dims(i_ax)}] = ismember(values{i_ax}, obj.(axes{i_ax}));
                if ops{dims(i_ax)} == 0
                    ops{dims(i_ax)} = [];
                end
            end

            obj = obj(ops{:});
            if nargout == 2
                indices = ops;
            end
        end

        function [obj, indices] = pickt(obj, axes, ranges)
            % Index by equality with tolerance
            % xarr.range(axis = [center tol])
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
            if nargout == 2
                indices = ops;
            end
        end

        function [obj, indices] = range(obj, axes, ranges)
            % Index by range membership, according to min(rng) <= coord & coord <= max(rng)
            % xarr.range(axis = [lower, upper])
            %   [upper, lower] also works, and either can be +/- Inf
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
            if nargout == 2
                indices = ops;
            end
        end


        function [obj, indices] = index(obj, axes, indices)
            % Index by logical or ordinal index
            % xarr.index(axis = [indices...])
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
            if nargout == 2
                indices = ops;
            end
        end

        function obj = rename(obj, old, new)
            % Rename xarray axis
            % xarr = xarr.rename(old, new)
            %   old     number or name of old axis
            %   new     name of new axis
            arguments
                obj xarray;
                old (1,1) {mustBeA(old, ["numeric", "string"])};
                new (1,1) string;
            end

            old = obj.convertdims(old);
            obj.axes(old) = new;
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
                mex = MException("xarray:axes", ...
                    "Unrecognized axis or property name(s) '%s'", names(~any(matches, 1)));
                throwAsCaller(mex);
            end
        end

        function iseq = dimequal(obj, other, dims)
            dims_a = obj.convertdims(dims);
            dims_b = other.convertdims(dims);
            iseq = false(1, length(dims));

            for i_dim = 1:length(dims)
                dim_a = dims_a(i_dim);
                dim_b = dims_b(i_dim);
                iseq(i_dim) = obj.axes(dim_a) == other.axes(dim_b) && ...
                    isequal(obj.coordinates{dim_a}, other.coordinates{dim_b});
            end
        end

        function indexcheck(obj, op)
            if length(op.Indices) ~= naxes(obj)
                mex = MException("xarray:indexing", ...
                    "Expected %d subscripts but got %d", ...
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

            new_coords = cell(1, length(obj.coordinates));
            for i_axis = 1:length(obj.coordinates)
                new_coords{i_axis} = obj.coordinates{i_axis}(op.Indices{i_axis});
            end
            obj.coordinates = new_coords;

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
                    error("Unable to assign sub-array: mismatch along %s", ...
                        obj.names(i_axis));
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
            old_sz = size(obj.data);
            obj.data(op.Indices{:}) = data;
            if ~isequal(size(obj.data), old_sz)
                error("Brace assignment changed xarray size from %d to %d", ...
                    old_sz, size(obj.data));
            end
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

        function obj = dotAssign(obj, operations, input)
            op = operations(1);
            axis = op.Name;
            dim = find(axis == obj.axes, 1);

            if isempty(dim) % new axis
                if ~isscalar(operations)
                    error("Invalid chained indexing while creating new axis");
                end
                if axis == ""
                    error("New axis must have non-empty name");
                end
                dim = length(obj.axes) + 1;
                obj.axes(dim) = axis;
                obj.coordinates{dim} = input;
            else % assign new values to existing axis
                if isscalar(operations)
                    obj.coordinates{dim} = input;
                elseif length(operations) == 2
                    obj.coordinates{dim}.(operations(2)) = input;
                end
            end
            
            % Validate newly assigned coordinates
            if ~isvector(obj.coordinates{dim})
                error("Coordinate variables must be vectors");
            end

            if length(obj.coordinates{dim}) ~= size(obj.data, dim)
                error("Mismatched axes: size(xarr, %d) = %d but coordinate has length %d", ...
                    dim, size(obj.data, dim), obj.coordinates{dim});
            end

            % Silently enforce column-vector coordinate variables by transposing row
            if size(obj.coordinates{dim}, 2) > 1
                obj.coordinates{dim} = obj.coordinates{dim}';
            end
        end

        function n = dotListLength(~, ~, ~)
            n = 1;
        end
    end

%% OBJECT DISPLAY OVERRIDE
    methods (Access = protected)
        function header = getHeader(obj)
            dimstr = matlab.mixin.CustomDisplay.convertDimensionsToString(obj);
            header = sprintf('%s %s with axes:', dimstr, class(obj));
        end

        function group = getPropertyGroups(obj)
            prop_list = [cellstr(obj.axes); 
                cellfun(@(c) c', obj.coordinates, UniformOutput = false)];
            prop_list = struct(prop_list{:});
            group = matlab.mixin.util.PropertyGroup(prop_list);
        end
    end
%
% %% MATH
%     % Internal templates for different "kinds" of operation
    methods (Static, Access = protected)
        % override an arithmetic operation based on common rules
        function obj = arithmetic(func, first, second)
            first_isx = isa(first, "xarray");
            second_isx = isa(second, "xarray");
            
            if first_isx && second_isx 
                % operate on two xarrays
                [order, common] = matchdims(first.axes, second.axes);

                iseq = dimequal(first, second, common);
                if any(~iseq)
                    mex = MException("xarray:mismatch", ...
                        "Common axes %s are not equal", ...
                        mat2str(common_names(~iseq)));
                    throwAsCaller(mex);
                end

                if isscalar(order)
                    order = [1 2];
                end
                second_data = double(second);
                second_data = permute(second_data, order);

                obj = first;
                obj.data = func(obj.data, second_data);
                obj.axes = [obj.axes second.axes(only_in_snd)];
                obj.coordinates = [obj.coordinates second.coordinates(only_in_snd)];
            else
                % operate on xarray and double (in either order)
                if first_isx
                    obj = first;
                    old_sz = size(obj.data);
                    obj.data = func(obj.data, second);
                elseif second_isx
                    obj = second;
                    old_sz = size(obj.data);
                    obj.data = func(first, obj.data);
                end

                if ~isequal(size(obj.data), old_sz)
                    mex = MException("xarray:arithmetic", ...
                        "The arithmetic operation expanded the object from %s to %s." + ...
                        " To create new dimensions by broadcasting element-wise operations," + ...
                        " Create an <xarray> instance to define the new axes.", ...
                        mat2str(old_sz), mat2str(size(obj.data)));
                    throwAsCaller(mex);
                end
            end
        end

        function [order, common] = matchdims(first, second)
            % [order, common] = matchdims(first, second)
            % first, second (string)    axis names from first and second arguments
            % axes common to both are tested for exact equality and ordered as in the first
            % axes only in the first are ordered as in the first
            % axes only in the second are appended to the end
            arguments
                first (1, :) string;
                second (1, :) string;
            end 

            matches = first' == second; % relies on xarray.axes forced to be a row
            % find axes in common
            [common_from_fst, common_from_snd] = find(matches);
            common = first.axes(common_from_fst);

            only_in_fst = find(~any(matches, 2));
            only_in_snd = find(~any(matches, 1));
            dims_total = size(matches, 1) + size(matches, 2) - length(common);
            new_dims = (size(matches, 2)+1):dims_total;

            % create permutation vector
            order = NaN(1, size(matches, 1));
            order(common_from_fst) = common_from_snd;
            order(only_in_fst) = new_dims; %#ok
            order = [order only_in_snd];
            
            assert(all(isfinite(permutation)), "Some permutation dimensions not assigned");
        end
    end
%
    % Overrides
    methods (Access = public)
        function obj = plus(a, b)
            obj = xarray.arithmetic(@plus, a, b);
        end

        function obj = minus(a, b)
            obj = xarray.arithmetic(@minus, a, b);
        end

        function obj = uminus(obj)
            obj = xarray.arithmetic(@minus, 0, obj);
        end

        function obj = times(a, b)
            obj = xarray.arithmetic(@times, a, b);
        end

        function obj = rdivide(a, b)
            obj = xarray.arithmetic(@rdivide, a, b);
        end

        function obj = power(a, b)
            obj = xarray.arithmetic(@power, a, b);
        end

        function obj = mtimes(a, b)
            obj = xarray.arithmetic(@times, a, b);
        end

        function obj = mrdivide(a, b)
            obj = xarray.arithmetic(@rdivide, a, b);
        end

        function obj = mpower(a, b)
            obj = xarray.arithmetic(@power, a, b);
        end

        function obj = sort(obj, dims, varargin)
            dims = obj.convertdims(dims);

            ops = repmat({':'}, 1, naxes(obj));
            for i_dim = 1:length(dims)
                dim = dims(i_dim);
                [~, ops{dim}] = sort(obj.coordinates{dim}, varargin{:});
            end
            obj = obj(ops{:});
        end

        function result = interp(obj, method, axes, values)
            % Interpolate over <xarray> dimensions
            arguments
                obj xarray;
                method (1,1) string;
            end

            arguments (Repeating)
                axes (1,1) string;
                values (1,:) double;
            end
            
            axes = string(axes);
            [interpolant, obj] = create_data_interpolant(obj, axes, method);

            result = interpolant(values{:});
            new_axes = obj.axes;
            new_coords = obj.coordinates;
            new_coords(1:length(values)) = values;
            result = xarray(result, new_axes, new_coords);
        end

        function interpfcn_h = interpolant(obj, axes, method, extrap)
            arguments
                obj xarray;
                axes (1,:) string;
                method (1,1) string = "linear";
                extrap (1,1) string = method;
            end

            [terp, obj] = create_data_interpolant(obj, axes, method, extrap);
            trailing_dims = (length(axes)+1):naxes(obj);
            trailing_coords = obj.coordinates(trailing_dims);
            % trailing_axes = obj.axes(trailing_dims);
            %
            % template_size = [ones(1, length(axes)), size(obj, trailing_dims)];
            % template_axes = []
            % template_coords = [repmat({1}, 1, length(axes)), obj.coordinates(trailing_dims)];
            %
            % lens = cellfun(@length, trailing_coords);
            % template = xarray(NaN(lens), trailing_axes, trailing_coords);

            interpfcn_h = @interpfcn;
            function result = interpfcn(varargin)
                if iscell(varargin{1})
                    gridvectors = varargin{1};
                else
                    gridvectors = varargin;
                end
                if length(gridvectors) ~= length(axes)
                    error("Expected %d grid vectors, %d specified", length(axes), length(gridvectors));
                end
                result = terp(gridvectors);
                result = xarray(result, obj.axes, [gridvectors trailing_coords]);
            end
        end

        function [interpolant, obj] = create_data_interpolant(obj, axes, method, extrap)
            % Create data interpolant
            % [interpolant, obj] = create_data_interpolant(obj, method, extrap, axes)
            %       Performs input checks and sorts specified axes
            arguments
                obj xarray;
                axes (1,:) string;
                method (1,1) string;
                extrap (1,1) string;
            end
            
            interp_dims = obj.convertdims(axes);

            err_id = "xarray:interpolate";
            lengths = cellfun(@length, obj.coordinates(interp_dims));
            scalars = axes(lengths == 1);
            if ~isempty(scalars)
                throwAsCaller(MException(err_id, "Dimension(s) %s are scalar." + ...
                    " Interpolation is not supported along scalar dimensions.", ...
                    mat2str(scalars)));
            end

            isnum = cellfun(@isnumeric, obj.coordinates(interp_dims));
            nonnum = axes(~isnum);
            if ~isempty(nonnum)
                throwAsCaller(MException(err_id, "Dimension(s) %s are non-numeric." + ...
                    " Interpolation is not supported along non-numeric dimensions.", ...
                    mat2str(nonnum)));
            end

            obj = obj.sort(interp_dims, "ascend");
            obj = obj.leading(interp_dims);
            interp_vectors = obj.coordinates(1:length(axes));
            interpolant = griddedInterpolant(interp_vectors, obj.data, method, extrap);
        end

        % function fun = interpolant(obj, dims, varargin)
        %     terp_dims = obj.convertdims(dims);
        %     terp_names = obj.axes(terp_dims);
        %     preserve_dims = 1:naxes(obj);
        %     preserve_dims(terp_dims) = [];
        %
        %     gridvectors = obj.coordinates(terp_dims);
        %     obj = sort(obj, terp_dims, "ascend");
        %     obj = permute(obj, [terp_dims preserve_dims]);
        %
        %     types = repmat("", 1, length(gridvectors));
        %     for i = 1:length(gridvectors)
        %         if isnumeric(gridvectors{i})
        %             types(i) = "numeric";
        %         elseif isduration(gridvectors{i})
        %             types(i) = "duration";
        %             gridvectors{i} = seconds(gridvectors{i});
        %         elseif isdatetime(gridvectors{i})
        %             types(i) = "datetime";
        %             gridvectors{i} = datenum(gridvectors{i});
        %         else
        %             error("Interpolation along non-numeric or time axes not supported");
        %         end
        %     end
        %
        %     terp = griddedInterpolant(gridvectors, double(obj), varargin{:});
        %     fun = @interpolate_and_assign;
        %
        %     function out = interpolate_and_assign(values)
        %         if length(values) ~= length(gridvectors)
        %         end
        %         values = varargin{2:2:end};
        %         if length(names) ~= length(values)
        %             error("Interpolant requires equal number of names and value vectors");
        %         end
        %         dims = obj.convertdims(names);
        %
        %         gridvalues = gridvectors;
        %         for i = 1:length(values)
        %             if isnumeric(values{i})
        %                 gridvalues{dims(i)} = values{i};
        %             elseif isduration(values{i});
        %                 gridvalues{dims(i)} = seconds(values{i});
        %             elseif isdatetime(values{i});
        %                 gridvalues{dims(i)} = datenum(values{i});
        %             else
        %                 error("Interpolation along non-numeric or time axes not supported");
        %             end
        %         end
        %         
        %     end
        % end
    end
end

