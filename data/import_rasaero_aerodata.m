%% Import aerodynamic data from RasAero 
% data = import_rasaero_table(path)
% Inputs
%   path    (string) path to .csv output from RasAero II
% Outputs
%   data    (xarray) xarray defining drag data, with axes ["mach", "aoa", "field"]
%                    "field" is the column name in the .csv table output

function aerodata = import_rasaero_aerodata(path)
    arguments (Input)
        path (1,1) string;
    end
    arguments (Output)
        aerodata (:, :, :) xarray {mustHaveAxes(aerodata, ["mach", "aoa", "field"])};
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


