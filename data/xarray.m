%% Array with named coordinate variables
% Barebones MATLAB implementation of python-xarray
% 
% EXAMPLES:
% % Create arrays
% xarr = xarray(randi(100, 10, 50, 3), ...
%     time = seconds(1:10), x = linspace(0, 10, 50), ...
%     name = ["first", "second", "third"])
% another_xarr = xarray(randi(100, 100, 50, 3), ...
%     time = seconds(11:110), x = linspace(0, 10, 50), ...
%     name = ["first", "second", "third"]) 
% % Join arrays by axis name
% xarr = cat("time", xarr, another_xarr)
% % Permute by name
% xarr = permute(xarr, ["name", "x", "time"])
% Select values by membership in range
% xarr = xarr.range(x = [2 5]) 
% % Select values by tolerance around value
% xarr = xarr.pickt(time = [seconds(4), seconds(1)], x = [4 0.5]) % select by tolerance around value
% Select subset by equality with specific value
% % Chain selections
% xarr = xarr.pick(name = "second") 
% % Chain selections
% another_xarr = another_xarr.range(x = [2 5]).pickt(time = [seconds(20) seconds(2)]) % chain selections

% NOTE FOR DEVELOPERS
% - orientation of <axes> and members of <coordinates> are deliberatly set to
% use matlab's broadcasting rules for certain operations

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
            % Length of coordinate variables must match size(data, n). Trailing
            % singleton dimensions are permitted (intended for later
            % concatenation)
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

            if length(axes) ~= length(coordinates)
                error("Specified %d axes but %d coordinates", ...
                    length(axes), length(coordinates));
            end

            obj.data = data;
            obj.axes = string(axes);
            obj.coordinates = coordinates;
        end

        function obj = set.data(obj, data)
            data = squeeze(data);
            if isvector(data) && size(data, 2) > 1
                data = data';
            end
            obj.data = data;
        end

        function obj = set.axes(obj, axes)
            err_id = "xarray:invalidAxes";
            if length(axes) ~= ndims(obj)
                throwAsCaller(MException( err_id, ...
                    "%d axes required, %d specified", ...
                    ndims(obj), length(axes)));
            end

            obj.axes = axes;
        end

        function obj = set.coordinates(obj, coords)
            err_id = "xarray:invalidAxes";
            if ndims(obj) ~= length(coords)
                throwAsCaller(MException( err_id, ...
                    "%d coordinates required, %d specified", ...
                    ndims(obj), length(coords) ));
            end
            for i = 1:length(coords)
                if ~isvector(coords{i})
                    throwAsCaller(MException( err_id, ...
                        "Coordinate %s is not a vector", i ));
                end
                if length(coords{i}) ~= size(obj, i)
                    throwAsCaller(MException( err_id, ...
                        "Coordinate %d are length %d, but data has size %d", ...
                        i, length(coords{i}), size(obj, i) ));
                end
                % always transpose row to column
                if size(coords{i}, 2) > 1
                    coords{i} = coords{i}';
                end
            end

            obj.coordinates = coords;
        end

        function dbl = double(obj)
            dbl = obj.data;
        end

        function nd = ndims(obj)
            nd = ndims(obj.data);
            if size(obj.data, 2) == 1
                nd = 1;
            end
        end

        function varargout = size(obj, varargin)
            if isempty(varargin)
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
                if dim > ndims(template)
                    error("Concatenation dimension %d too large (%d-D)." + ...
                        " To create a new dimension, specify the axis as a string scalar")
                end
            else
                name = dim;
                dim = find(template.axes == name);

                if isempty(dim)
                    dim = ndims(template) + 1;
                end
            end
            
            to_compare = 1:ndims(template);

            % create new coordinate axis
            if dim <= ndims(template) % existing dimension
                coords = cellfun(@(arry) arry.coordinates{dim}, ...
                    arrays, UniformOutput = false);
                coords = vertcat(coords{:});
                to_compare(dim) = [];
            else % new dimension defaults to 1:n
                coords = 1:length(arrays);
            end


            for i_obj = 2:length(arrays)
                if any(~dimequal(template, arrays{i_obj}, to_compare))
                    error("Mismatched dimensions at position %d", i_obj);
                end
            end

            datas = cellfun(@(arry) arry.data, arrays, UniformOutput = false);
            template.data = cat(dim, datas{:});
            template.coordinates{dim} = coords;
            obj = template;
        end


        function obj = permute(obj, order)
            % Re-order axes by name or number
            % xarr = permute(xarr, [dim1 dim2 ...]);

            dimorder = obj.convertdims(order);

            if length(dimorder) ~= ndims(obj)
                error("Dimension order must have exactly one entry per axis:\n%s", ...
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

        function obj = pick(obj, axes, values)
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
            ops = repmat({':'}, 1, ndims(obj));

            for i_ax = 1:length(axes)
                [~, ops{dims(i_ax)}] = ismember(values{i_ax}, obj.(axes{i_ax}));
                if ops{dims(i_ax)} == 0
                    ops{dims(i_ax)} = [];
                end
            end

            obj = obj(ops{:});
        end

        function obj = pickt(obj, axes, ranges)
            % Index by equality with tolerance
            % xarr.range(axis = [center tol])
            arguments
                obj xarray
            end
            arguments (Repeating)
                axes (1, 1) string;
                ranges (1, 2);
            end
            
            ops = repmat({':'}, 1, ndims(obj));
            dims = obj.convertdims(axes{:});

            for i_ax = 1:length(axes)
                coord = obj.(axes{i_ax});
                rang = ranges{i_ax}(1) + ranges{i_ax}(2) * [-1 1];
                ops{dims(i_ax)} = rang(1) <= coord & coord <= rang(2);
            end

            obj = obj(ops{:});
        end

        function obj = range(obj, axes, ranges)
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

            ops = repmat({':'}, 1, ndims(obj));
            dims = obj.convertdims(axes{:});

            for i_ax = 1:length(axes)
                coord = obj.(axes{i_ax});
                ops{dims(i_ax)} = min(ranges{i_ax}) <= coord & coord <= max(ranges{i_ax});
            end

            obj = obj(ops{:});
        end


        function obj = index(obj, axes, indices)
            % Index by logical or ordinal index
            % xarr.index(axis = [indices...])
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
                ops{dims(i_ax)} = indices{i_ax};
            end

            obj = obj(ops{:});
        end

        function obj = rename(obj, old, new)
            arguments
                obj xarray;
                old (1,1) {mustBeA(old, ["numeric", "string"])};
                new (1,1) string;
            end
            % Rename xarray axis
            % xarr = xarr.rename(old, new)
            %   old     number or name of old axis
            %   new     name of new axis
            old = obj.convertdims(old);
            obj.axes(old) = new;
        end
    end


    methods (Access = public)
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
            if length(op.Indices) ~= ndims(obj)
                mex = MException("xarray:indexing", ...
                    "Expected %d subscripts but got %d", ...
                    ndims(obj), length(op.Indices));
                throwAsCaller(mex);
            end
        end

%% INDEXING OVERRIDES
        function obj = parenReference(obj, operations)
            % Input processing and checking
            op = operations(1);
            indexcheck(obj, op);

            % Index into data using all indices, index into each axis using its index
            new_data = obj.data.(op);
            sizes = size(new_data);
            nonscalar = sizes ~= 1;
            obj.data = new_data;
            obj.axes = obj.axes(nonscalar);

            new_coords = cell(1, length(obj.coordinates));
            for i_axis = 1:length(obj.coordinates)
                new_coords{i_axis} = obj.coordinates{i_axis}(op.Indices{i_axis});
            end
            obj.coordinates = new_coords(nonscalar);

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
            for i_axis = 1:ndims(obj)
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

            for i_axis = 1:ndims(obj)
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

        function obj = dotAssign(obj, operations, input)
            op = operations(1);
            axis = op.Name;
            axis_idx = convertdims(obj, axis);

            if isscalar(operations)
                obj.coordinates{axis_idx} = input;
            elseif length(operations) == 2
                obj.coordinates{axis_idx}.(operations(2)) = input;
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
                % obj.coordinates];
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
                % axes common to both are tested for exact equality and ordered as in the first
                % axes only in the first are ordered as in the first
                % axes only in the second are appended to the end

                matches = first.axes' == second.axes; % relies on xarray.axes forced to be a row
                % find axes in common
                [common_from_fst, common_from_snd] = find(matches);
                common_names = first.axes(common_from_fst);

                iseq = dimequal(first, second, common_names);
                if any(~iseq)
                    mex = MException("xarray:mismatch", ...
                        "Common axes %s are not equal", ...
                        mat2str(common_names(~iseq)));
                    throwAsCaller(mex);
                end

                % create permutation vector
                only_in_fst = find(~any(matches, 2));
                only_in_snd = find(~any(matches, 1));
                dims_total = size(matches, 1) + size(matches, 2) - length(common_names);
                new_dims = (size(matches, 2)+1):dims_total;

                permutation = NaN(1, size(matches, 1));
                permutation(common_from_fst) = common_from_snd;
                permutation(only_in_fst) = new_dims; %#ok
                permutation = [permutation only_in_snd];

                assert(all(isfinite(permutation)), "Some permutation dimensions not assigned");

                if isscalar(permutation)
                    permutation = [1 2];
                end
                second_data = double(second);
                second_data = permute(second_data, permutation);

                obj = first;
                obj.data = func(obj.data, second_data);
                obj.axes = [obj.axes second.axes(only_in_snd)];
                obj.coordinates = [obj.coordinates second.coordinates(only_in_snd)];
            else
                % operate on xarray and double (in either order)
                if first_isx
                    obj = first;
                elseif second_isx
                    obj = second;
                end

                old_nd = ndims(obj);
                obj.data = func(double(first), double(second));
                
                if ndims(obj) > old_nd
                    mex = MException("xarray:arithmetic", ...
                        "The arithmetic operation expanded the object to %d dimensions." + ...
                        " To create new dimensions by broadcasting element-wise operations," + ...
                        " Create an <xarray> instance to define the new axes.", ndims(obj));
                    throwAsCaller(mex);
                end
            end
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

            ops = repmat({':'}, 1, ndims(obj));
            for i_dim = 1:length(dims)
                dim = dims(i_dim);
                [~, ops{dim}] = sort(obj.coordinates{dim}, varargin{:});
            end
            obj = obj(ops{:});
        end

        % function fun = interpolant(obj, dims, varargin)
        %     terp_dims = obj.convertdims(dims);
        %     terp_names = obj.axes(terp_dims);
        %     preserve_dims = 1:ndims(obj);
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

