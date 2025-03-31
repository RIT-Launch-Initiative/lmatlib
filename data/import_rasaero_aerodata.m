%% lookups = import_rasaero_table(path, format = "grid")
% Import aerodynamic data from RasAero
% INPUTS
%   path:   path to .csv, .xlsx, or other table file
%   format: (optional) format of output structure (default: grid)
% OUTPUTS
%   lookups: structure with following fields /
%       Ma:     mach number values
%       alpha:  angle of attack values
%   + One field for each additonal variable in the table (CD, CN, CA, ...):
%       For "grid" format, this is an length(Ma) by length(alpha) matrix with
%       this variable's lookup values.
%       For "terp" format, this is a griddedInterpolant object that can be
%       queried with CD(Ma, alpha).

function aerodata = import_rasaero_aerodata(path)
    arguments
        path (1,1) string;
    end
    data = readtable(path, VariableNamingRule = "preserve");
    names = string(data.Properties.VariableNames);
    [Ma, alpha, schema] = unstack_indices(data.Mach, data.Alpha);
    data_names = names(~contains(names, ["Mach", "Alpha"]));
    aerodata = xarray(zeros(length(Ma), length(alpha), length(data_names)), ...
        ["mach", "aoa", "field"], {Ma, deg2rad(alpha), data_names});

    for var = data_names
        aerodata.pick{"field", var} = data{:, var}(schema);
    end
end

% use unstack() to find 
% returns 
%   x_1, x_2: The unique values of x_1 and x_2
%   imat: the destination for each row when unstacked into length(x_1) by length(x_2) matrix
function [x_1, x_2, imat] = unstack_indices(x_1, x_2)
    t = table(x_1, x_2);
    t.Indices = (1:height(t))';
    t = unstack(t, "Indices", "x_2", ...
        VariableNamingRule = "preserve");
    t = fillmissing(t, "nearest");

    x_1 = t.x_1;
    names = t.Properties.VariableNames;
    imat = t{:, names ~= "x_1"};
    x_2 = str2double(names(names ~= "x_1"));
end


