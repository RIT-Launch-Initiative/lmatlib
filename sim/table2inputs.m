%% Create simulation input array from table of variables
% simins = table2inputs(simin, table)
% Inputs
%   simin   Simulink.SimulationInput object
%   table   Table of variables to assign
% Outputs
%   simins  Simulink.SimulationInput object array
%           Each table row produces one object

function simins = table2inputs(simin, table)
    n = height(table);
    c = width(table);
    names = string(table.Properties.VariableNames)';
    simins(1:n) = simin;
    
    for row = 1:n
        for col = 1:c
            simins(row) = simins(row).setVariable(names(c), table{row, col});
        end
    end
end
