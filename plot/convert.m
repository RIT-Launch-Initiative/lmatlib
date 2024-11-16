% Convert all occurences of one unit in a timetable to another unit, using the supplied function
% converted_table = convert(tt, from, to, fun)
%   tt:     input timetable
%   from:   unit to convert
%   to:     unit to convert all "from" to
%   fun:    function to use to convert
%           fun output must be the same size as its input
%
% Example: convert(table, "deg", "rad", @deg2rad)
% Example: convert(table, "ft/s", "km/s", @(col) col ./ 3281)
% ! NOTE: Operates by string matching - will not interpret "deg" as the same unit as the degree symbol

function tt = convert(tt, from, to, fun)
    arguments 
        tt timetable;
        from (1,1) string;
        to (1,1) string;
        fun function_handle;
    end

    from_cols = tt.Properties.VariableUnits == from;
    
    if ~any(from_cols)
        error("Unit ""%s"" does not exist in the input timetable", from);
    end

    tt{:, from_cols} = fun(tt{:, from_cols});
    tt.Properties.VariableUnits(from_cols) = to;
end
