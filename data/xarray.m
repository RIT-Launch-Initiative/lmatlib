%% Array with named coordinate variables
% 
% Some data are highly dimensional. Weather data, for example, can have as many
% as five: (time, lat, lon, height, variable). Using this data as a raw 5-D
% double requires separately tracking the order of each dimension and the
% coordinate variable associated with each dimension. <xarray> associates named
% coordinate variables to the dimensions of a double so that code operating on
% the array more directly expresses to the author's intent. 
%
% This class is inspired by, but is not written or endorsed by the authors of
% the original Python xarray.
% 
% xarray(data, axis1 = values1, axis2 = values2, ...) creates an xarray. 
%   The resulting xarray will have dimensions named (axis1, axis2) and each
%   dimension will have coordinate variables (values1, values2, ...).
% Each set of values must be a row or column vector. All coordinates are
% implicitly converted to column-vectors. Taken in order, they must have the
% same size as the data -- 
% length(values1) == size(data, 1); 
% length(values2) == size(data, 2);
% 
% SELECTION
% Currently, four accessors are supported: "pick", "pickt", "range", and
% "index". They select elements of <xarray> by finding an axis by name instead
% of index.
% xarr.pick(axis1 = [value, another_value], ...) 
%   returns a slice of <xarr> where the value <axis1> matches the selected
%   values---in the order specified. 
% xarr.pickt(axis1 = [center tol], ...) 
%   returns a slice of <xarr> where the value of <axis1> is within <tol> of
%   <center>.
% xarr.range(axis1 = [start end], ...)
%   returns a slice of <xarr> where the value of <axis1> is 
% xarr.range(axis1 = indices, ...)
%   returns a slice of <xarr> by index
% 
% In general, brace-indexing applies to the underlying data (in the same way it
% does for tables), but Name=Value syntax only works for paren-operations in MATLAB.
% out = xarr.range(axis1 = indices, ...)
%   becomes
% out = xarr.range{"axis1", indices, ...}
%   and performs the same indexing operation, but returns the raw data.
% 
% ASSIGNMENT
% xarr.pick(___) = other_xarr
% xarr.pick{___} = data
%   assigns the data of <input> to the slice that would be returned by the
%   corresponding accessor.
% 
% ARITHMETIC
% xarray only supports element-wise operations. 
% <xarray> +/* <double> 
%   applies the operation element-wise to the data. However, this may not
%   create or change dimensions: if the result is a different size from the
%   original value of data, this is an error.
% <xarray> +/* <xarray>
%   aligns common dimensions of each <xarray>, tests them for exact equality,
%   and broadcasts operations across non-shared dimensions. 
% <xarray: x, y, z> * <xarray: t, x, z> = <xarray: x, y, z, t>, 
%   if the values of the x and z axes x and z are exactly equal in both.
% 
% SUMMARIZING FUNCTIONS
% mx = mean(xarray, "x")
% The other summarizing funcitons (min, max, vecnorm, ...) take the same
% additional arguments as their numeric equivalents, but the additional
% arguments are required.
% 

% NOTES FOR DEVELOPERS
% Paren-indexing applies corresponding indices to the axes and coordinates for
% reference and assignment. This forms the base for access operations.
% out = xarr(1:2, :)
%   returns an xarray with the first two values along the first axis and all
%   values along the second axis.
% xarr(1:2, :) = input     
%   compares the dimensions of xarr(1:2, :) and input, assigning the input data
%   to xarr if and only if the dimensions exactly match.
%
% Linear indexing is not supported.
%
% Axis and data access are both implemented by overriding dotReference/Assign,
% because that was the simplest way to support assignment using the custom
% accessors - if they were methods at the object level, assignment to their
% return value would not be valid MATLAB. The dot operations therefore have two
% branches: one for axis access/assignment, one for data access/assignment.
% Axis access finds the axis of the appropriate name. Data access matches the
% specified axis names to numeric dimensions and uses paren/brace
% Reference/Assign. 

classdef xarray < matlab.mixin.indexing.RedefinesDot ...
        & matlab.mixin.indexing.RedefinesParen ...
        & matlab.mixin.indexing.RedefinesBrace ...
        & matlab.mixin.CustomDisplay

    properties (GetAccess = public, SetAccess = protected)
        axes (1, :) string = [];
        coordinates (1, :) cell = {};
        data = [];
    end

    properties (Access = protected)
        access_labels (1, :) string {mustBeNonempty}= ...
            ["pick", "pickt", "range", "index"];
        access_fcns (1, :) cell {mustBeNonempty}= ...
            {@acs_pick, @acs_pickt, @acs_range, @acs_index};
    end

    methods
        function obj = xarray(data, axes, coordinates)
            % Construct array with named coordinate variables
            % xarr = xarray(data, axis1 = coordinates1, axis2 = coordinates2[, ...])
            %   data                double array
            %   axis1 ...           string axis name
            %   coordinates1 ...    coordinate variables (vectors, any type)
            %
            % Coordinate variables must
            %   - Be vectors
            %   - Match the size of the data in that dimension
            arguments
                data {mustBeNumericOrLogical} = [];
            end

            arguments (Repeating)
                axes (1, :) string;
                coordinates (1, :);
            end

            if nargin == 0
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

        function [tf, dims] = ismatrix(da)
            dims = find(size(da) > 1);
            tf = length(dims) == 2;
        end

        function [tf, dims] = isvector(da)
            dims = find(size(da) > 1);
            tf = isscalar(dims);
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

        function obj = align(obj, order)
            % xarr = xarr.align(order) 
            leading_order = obj.convertdims(order);
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

        function mustHaveAxes(xarr, axes, types)
            % xarray argument validation function
            % checks that the xarray input has appropriate axis names (optionally, types) for the code
            % mustHaveAxes(xarr, axes[, types])
            %   Inputs
            %   xarr    (xarray)    xarray to validate
            %   axes    (string)    axes to check for
            %   types   (string)    (optional) types to check for (through <isa>)
            arguments
                xarr xarray;
                axes (1, :) string {mustBeNonempty};
                types (1, :) string = [];
            end
            
            err_id = "xarray:invalidAxes";
            % name = inputname(1);
            notfound = setdiff(axes, xarr.axes);
            if ~isempty(notfound)
                mex = MException(err_id, "Invalid xarray input: required axe(s) %s not present", ...
                    mat2str(notfound));
                throwAsCaller(mex);
            end
            
            if isempty(types)
                return;
            elseif length(types) ~= length(axes)
                error("Incorrect use of validator: %d axe(s) but %d type(s)", ...
                    length(axes), length(types))
            else
                correct = false(1, length(axes));
                for i_axis = 1:length(axes)
                    correct(i_axis) = isa(xarr.(axes(i_axis)), types(i_axis));
                end
                if any(~correct)
                    mex = MException(err_id, "Invalid xarray input: axe(s) %s must be %s", ...
                        mat2str(axes(~correct)), mat2str(types(~correct)));
                    throwAsCaller(mex);
                end
                
            end
        end

        function dims = convertdims(obj, varargin)
            % Convert named or numbered axes to numeric dimensions
            % dims = convertdims(xarr, [axis1 axis2])
            % dims = convertdims(xarr, axis1, axis2)
            arguments (Output)
                dims (1,:) double {mustBeInteger, mustBePositive};
            end
            
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

        function gh = plot(obj, varargin)
            arguments
                obj xarray;
            end
            arguments (Repeating)
                varargin
            end

            [isv, dim] = isvector(obj);
            assert(isv, "<xarray> input to <plot> must have exactly one non-scalar dimension");

            gh = plot(obj.coordinates{dim}, obj.data(:), varargin{:});
            if gh.Parent.XLabel.String == ""
                gh.Parent.XLabel.String = obj.axes(dim);
            end
        end

        function gh = imagesc(obj, opts)
            arguments
                obj xarray;
                opts.cmap (1,1) string = "parula";
                opts.clabel (1,1) string = missing;
            end

            [ism, dim] = ismatrix(obj);
            assert(ism, "<xarray> input to <imagesc> must have exactly two non-scalar dimensions");
            obj = obj.align(dim);
            gh = imagesc(obj.coordinates{2}, obj.coordinates{1}, squeeze(obj.data));

            ax = gh.Parent;
            ax.YDir = "normal";

            if ax.XLabel.String == ""
                ax.XLabel.String = obj.axes(2);
            end
            if ax.YLabel.String == ""
                ax.YLabel.String = obj.axes(1);
            end

            colormap(ax, opts.cmap);
            if ~ismissing(opts.clabel)
                cb = colorbar(ax);
                cb.Label.String = opts.clabel;
            end
            axis(ax, "tight");
        end

        function gh = surf(obj, opts)
            arguments
                obj xarray;
                opts.cmap (1,1) string = "parula";
                opts.clabel (1,1) string = missing;
            end

            [ism, dim] = ismatrix(obj);
            assert(ism, "<xarray> input to <imagesc> must have exactly two non-scalar dimensions");
            obj = obj.align(dim);
            gh = surf(obj.coordinates{2}, obj.coordinates{1}, squeeze(obj.data));

            ax = gh.Parent;

            if ax.XLabel.String == ""
                ax.XLabel.String = obj.axes(2);
            end
            if ax.YLabel.String == ""
                ax.YLabel.String = obj.axes(1);
            end

            colormap(ax, opts.cmap);
            if ~ismissing(opts.clabel)
                cb = colorbar(ax);
                cb.Label.String = opts.clabel;
            end
            axis(ax, "tight");
        end

        function gh = contour(obj, opts)
            arguments (Input)
                obj xarray
                opts.cmap (1,1) string = "parula";
                opts.clabel (1,1) string = missing;
                opts.levels (1,:) double = [];
            end

            [ism, dim] = ismatrix(obj);
            assert(ism, "<xarray> input to <contour> must have exactly two non-scalar dimensions");
            obj = obj.align(dim);
            if isempty(opts.levels)
                [~, gh] = contour(obj.coordinates{2}, obj.coordinates{1}, ...
                    squeeze(obj.data));
            else
                [~, gh] = contour(obj.coordinates{2}, obj.coordinates{1}, ...
                    squeeze(obj.data), opts.levels);
            end

            ax = gh.Parent;
            if ax.XLabel.String == ""
                ax.XLabel.String = obj.axes(2);
            end
            if ax.YLabel.String == ""
                ax.YLabel.String = obj.axes(1);
            end

            colormap(ax, opts.cmap);
            if ~ismissing(opts.clabel)
                cb = colorbar(ax);
                cb.Label.String = opts.clabel;
            end
            axis(ax, "tight");
        end

        function gh = contourf(obj, opts)
            arguments (Input)
                obj xarray
                opts.cmap (1,1) string = "parula";
                opts.clabel (1,1) string = missing;
                opts.levels (1,:) double = [];
            end

            [ism, dim] = ismatrix(obj);
            assert(ism, "<xarray> input to <contour> must have exactly two non-scalar dimensions");
            obj = obj.align(dim);
            if isempty(opts.levels)
                [~, gh] = contourf(obj.coordinates{2}, obj.coordinates{1}, ...
                    squeeze(obj.data));
            else
                [~, gh] = contourf(obj.coordinates{2}, obj.coordinates{1}, ...
                    squeeze(obj.data), opts.levels);
            end

            ax = gh.Parent;
            if ax.XLabel.String == ""
                ax.XLabel.String = obj.axes(2);
            end
            if ax.YLabel.String == ""
                ax.YLabel.String = obj.axes(1);
            end

            colormap(ax, opts.cmap);
            if ~ismissing(opts.clabel)
                cb = colorbar(ax);
                cb.Label.String = opts.clabel;
            end
            axis(ax, "tight");
        end
    end

    methods (Access = protected)
%% INDEXING HELPER FUNCTIONS
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

        function indices = acs_pick(obj, dim, param)
            coord = obj.coordinates{dim};
            [present, indices] = ismember(param, coord);
            if any(~present)
                throwAsCaller(MException("xarray:indexing", "Axis '%s' does not contain %s", ...
                    obj.axes(dim), mat2str(param(~present))));
            end
        end

        function indices = acs_pickt(obj, dim, param)
            coord = obj.coordinates{dim};
            indices = ((param(1) - param(2)) < coord) & (coord < (param(1) + param(2)));
        end

        function indices = acs_range(obj, dim, param)
            coord = obj.coordinates{dim};
            indices = (min(param) < coord) & (coord < max(param));
        end

        function indices = acs_index(~, ~, param)
            indices = param;
        end

        function indices = access(obj, fcn, axes, params)
            arguments
                obj xarray
                fcn (1,1) function_handle;
            end 
            arguments (Repeating)
                axes (1,1) string;
                params;
            end
            
            indices = repmat({':'}, 1, naxes(obj));
            dims = obj.convertdims(axes{:});
            for i_dim = 1:length(dims)
                dim = dims(i_dim);
                indices{dim} = fcn(obj, dim, params{i_dim});
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
            arguments
                obj xarray;
                operations;
                subobj xarray;
            end
            
            % Input processing and checking
            if ~isscalar(operations)
                error("Chained paren-assignment not supported")
            end
            op = operations(1);
            indexcheck(obj, op);
            dest = obj.(op);

            if ~isequal(size(dest), size(subobj))
                error("Unable to perform assignment because the left side is %s and the right side is %s", ...
                    string(size(dest)).join("-by"), string(size(subobj)).join("-by-"));
            end

            for i_axis = 1:naxes(dest)
                if dest.axes(i_axis) ~= subobj.axes(i_axis)
                    error("Mismatched names at dimension %d", i_axis);
                end
                if ~isequal(dest.coordinates{i_axis}, subobj.coordinates{i_axis})
                    error("Mismatched coordinate values at dimension %d", i_axis);
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
        function result = dotReference(obj, operations)
            op = operations(1);

            i_accessor = find(op.Name == obj.access_labels);
            if isempty(i_accessor) 
                % axis access
                dim = obj.convertdims(op.Name);
                result = obj.coordinates{dim};
                if ~isscalar(operations)
                    result = result.(operations(2:end));
                end

            else
                % indexor access
                if isscalar(operations) 
                    error("Accessors must have additional brace or paren-indices")    
                end
                indices = obj.access(obj.access_fcns{i_accessor}, operations(2).Indices{:});
                if operations(2).Type == "Paren"
                    result = obj(indices{:});
                elseif operations(2).Type == "Brace"                
                    result = obj{indices{:}};
                else
                    error("Accessors must have indices specified as name-value pairs")    
                end

                if length(operations) > 2
                    result = result.(operations(3:end));
                end
            end
        end

        function obj = dotAssign(obj, operations, input)
            op = operations(1);
            i_accessor = find(op.Name == obj.access_labels);
            if isempty(i_accessor)
                % axis creation or acccess
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
                        dim, size(obj.data, dim), length(obj.coordinates{dim}));
                end

                % Silently enforce column-vector coordinate variables by transposing row
                if size(obj.coordinates{dim}, 2) > 1
                    obj.coordinates{dim} = obj.coordinates{dim}';
                end
            else
                if length(operations) ~= 2 
                    error("Accessor-based assignment must be of the form\n\t<array>.<access>(<parameters>)")    
                end

                indices = obj.access(obj.access_fcns{i_accessor}, operations(2).Indices{:});
                if operations(2).Type == "Paren"
                    obj(indices{:}) = input;
                elseif operations(2).Type == "Brace"                
                    obj{indices{:}} = input;
                else
                    error("Accessors must have indices specified as name-value pairs")    
                end
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

%% MATH
    % Internal templates for different "kinds" of operation
    methods (Static, Access = protected)
        % override a summarizing operation
        function obj = summarize(func, obj, i_dims, args)
            err_id = "xarray:summarize";
            sz = size(obj);
            if length(args) < i_dims
                throwAsCaller(MException(err_id, ...
                    "Additional arguments to summarizing functions specifying dim/vecdim " + ...
                    "are required, in the same position as the base MATLAB equivalent. " + ...
                    "See\n\nhelp %s", func2str(func)));
            end
            if args{i_dims} == "all"
                args{i_dims} = 1:naxes(obj);
            else
                args{i_dims} = obj.convertdims(args{i_dims});
            end

            sz(args{i_dims}) = [];

            obj.data = func(obj.data, args{:});
            obj.data = reshape(obj.data, [sz 1]); % get rid of the singleton
            obj.axes(args{i_dims}) = [];
            obj.coordinates(args{i_dims}) = [];
        end

        % override an arithmetic operation based on common rules
        function obj = arithmetic(func, first, second)
            first_isx = isa(first, "xarray");
            second_isx = isa(second, "xarray");
            
            if first_isx && second_isx 
                % operate on two xarrays
                matches = first.axes' == second.axes; % relies on xarray.axes forced to be a row
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
                
                assert(all(isfinite(order)), "Some permutation dimensions not assigned");

                iseq = dimequal(first, second, common);
                if any(~iseq)
                    mex = MException("xarray:mismatch", ...
                        "Common axes %s are not equal", ...
                        mat2str(common(~iseq)));
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
    end

    % Math Overrides
    methods (Access = public)
%% Arithmetic
        % I don't have a better idea of how to do these unfortunately
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
        function obj = lt(a, b)
            obj = xarray.arithmetic(@lt, a, b);
        end
        function obj = gt(a, b)
            obj = xarray.arithmetic(@gt, a, b);
        end
        function obj = le(a, b)
            obj = xarray.arithmetic(@le, a, b);
        end
        function obj = ge(a, b)
            obj = xarray.arithmetic(@ge, a, b);
        end
        function obj = eq(a, b)
            obj = xarray.arithmetic(@eq, a, b);
        end
        function obj = ne(a, b)
            obj = xarray.arithmetic(@ne, a, b);
        end
        function obj = atan2(a, b)
            obj = xarray.arithmetic(@atan2, a, b);
        end

%% Summarizing functions
        function obj = mean(obj, varargin)
            obj = xarray.summarize(@mean, obj, 1, varargin);
        end
        function obj = median(obj, varargin)
            obj = xarray.summarize(@median, obj, 1, varargin);
        end

        function obj = var(obj, varargin)
            obj = xarray.summarize(@var, obj, 2, varargin);
        end
        function obj = std(obj, varargin)
            obj = xarray.summarize(@std, obj, 2, varargin);
        end
        function obj = vecnorm(obj, varargin)
            obj = xarray.summarize(@vecnorm, obj, 2, varargin);
        end
        function obj = min(obj, varargin)
            obj = xarray.summarize(@min, obj, 2, varargin);
        end
        function obj = max(obj, varargin)
            obj = xarray.summarize(@max, obj, 2, varargin);
        end
        function obj = sum(obj, varargin)
            obj = xarray.summarize(@sum, obj, 1, varargin);
        end
        function obj = rms(obj, varargin)
            obj = xarray.summarize(@rms, obj, 1, varargin);
        end

        function obj = sqrt(obj)
            obj.data = sqrt(obj.data);
        end
        function obj = abs(obj)
            obj.data = abs(obj.data);
        end
        function obj = log(obj)
            obj.data = log(obj.data);
        end
        function obj = sin(obj)
            obj.data = sin(obj.data);
        end
        function obj = cos(obj)
            obj.data = cos(obj.data);
        end
        function obj = tan(obj)
            obj.data = tan(obj.data);
        end
        function obj = asin(obj)
            obj.data = asin(obj.data);
        end
        function obj = acos(obj)
            obj.data = acos(obj.data);
        end
        function obj = atan(obj)
            obj.data = atan(obj.data);
        end

        function obj = isnan(obj)
            obj.data = isnan(obj.data);
        end
        function obj = isfinite(obj)
            obj.data = isfinite(obj.data);
        end
        function obj = any(obj, varargin)
            obj = xarray.summarize(@any, obj, 1, varargin);
        end
        function obj = all(obj, varargin)
            obj = xarray.summarize(@all, obj, 1, varargin);
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

        function obj = interp(obj, axes, values)
            arguments
                obj xarray;
            end
            arguments (Repeating)
                axes (1,1) string;
                values (1,:) double;
            end

            axes = [axes{:}];
            obj = obj.align(axes).sort(axes, "ascend");
            F = griddedInterpolant(obj.coordinates(1:length(axes)), obj.data, "linear", "nearest");
            result_data = F(values);
            new_size = size(result_data, length(axes)+1:ndims(result_data));
            obj.data = reshape(result_data, [new_size 1]);
            obj.axes(1:length(axes)) = [];
            obj.coordinates(1:length(axes)) = [];
        end
    end
end

