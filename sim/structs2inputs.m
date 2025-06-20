%% Assign structures' variables to simulation inputs
% simin = structs2inputs(simin, ...)
% Inputs
%   simin   Simulink.SimulationInput object -or- simulation name
%   ...     Any number of structures as separate arguments 
%           NOTE: If any of these structures are non-scalar, they must have the same
%           number of elements so the lengths stay compatible
%           NOTE: Scalar structures are assigned first, then non-scalar structures
%           NOTE: The usual use will be many scalar structures with different
%           parameter groups, and maybe one non-scalar array for sweeps or
%           Monte Carlo inputs
% Outputs
%   simin   Simulink.SimulationInput object(s) with assigned variables
%
% % Example
% rocket.A_ref = 23.1e-2;
% rocket.h_0 = 2300;
% rocket.v_0 = 300;
% simin = structs2inputs("sim_2dof", rocket); % the string "sim_2dof" is converted to a 
%                                           % SimulationInput and A_ref, h_0, v_0 appear 
%                                           in the model's workspace when it runs
% 
% config.dt = 0.01;
% simin = structs2inputs(simin, config); % can be used on an existing SimulationInput to add more variables
% simout = sim(simin);
% 
% Equivalently, in one line:
% simin = struct2inputs("sim_2dof", rocket, config) 

function simins = structs2inputs(simin, params)
    arguments (Input)
        simin {mustBeA(simin, ["string", "Simulink.SimulationInput"])};
    end
    arguments (Input, Repeating)
        params struct;
    end
    arguments (Output)
        simins Simulink.SimulationInput;
    end

    if isstring(simin)
        simin = Simulink.SimulationInput(simin);
    end

    lengths = cellfun(@numel, params);
    nonscalars = lengths(lengths ~= 1);
    if ~isscalar(nonscalars) && ~isempty(nonscalars)
        error("If there are multiple non-scalar structure inputs, " + ...
            "they must have the same number of elements");
    end

    i_scalars = find(lengths == 1);

    for i_struct = i_scalars
        this_params = params{i_struct};
        fields = string(fieldnames(this_params));
        for fn = 1:length(fields)
            name = fields(fn);
            simin = simin.setVariable(name, this_params.(name));
        end
    end

    if isempty(nonscalars)
        num_out = 1;
    else
        num_out = nonscalars;
    end

    simins(1:num_out) = simin;
    i_nonscalars = find(lengths ~= 1);

    for i_struct = i_nonscalars
        this_params = params{i_struct};
        fields = string(fieldnames(this_params));
        for i_sim = 1:numel(simins)
            for fn = 1:length(fields)
                name = fields(fn);
                simins(i_sim) = simins(i_sim).setVariable(name, ...
                    this_params(i_sim).(name));
            end
        end
    end
end
