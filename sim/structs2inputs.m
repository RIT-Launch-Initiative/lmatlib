%% Assign structures' variables to simulation inputs
% simin = structs2inputs(simin, params)
% Inputs
%   simin   Simulink.SimulationInput object -or- simulation name
%   params  Structure (array) with parameters
% Outputs
%   simin   Simulink.SimulationInput object(s) with assigned variables
%
% % Example
% rocket.A_ref = 23.1e-2;
% rocket.h_0 = 2300;
% rocket.v_0 = 300;
% simin = struct2input("sim_2dof", rocket); % the string "sim_2dof" is converted to a 
%                                           % SimulationInput and A_ref, h_0, v_0 appear 
%                                           in the model's workspace when it runs
% 
% config.dt = 0.01;
% simin = struct2input(simin, config); % can be used on an existing SimulationInput to add more variables
% simout = sim(simin);

function simins = structs2inputs(simin, params)
    arguments
        simin
        params struct
    end
    if isstring(simin)
        simin = Simulink.SimulationInput(simin);
    end

    simins(1:numel(params)) = simin;

    fields = string(fieldnames(params));
    for i_sim = 1:numel(simins)
        for fn = 1:length(fields)
            name = fields(fn);
            simins(i_sim) = simins(i_sim).setVariable(name, params(i_sim).(name));
        end
    end
end
