set(groot, "DefaultAxesNextPlot", "add");
clear; close all;
% NOTE - Aerospace Blockset required for the simulation file


%% Flight data import
SL = 1401; % [m] Ground level above sea level
openrocket = import_openrocket_csv("data/OMEN_OR_Output.csv");
openrocket.ASL = SL + openrocket.("Altitude");
openrocket.Properties.VariableUnits("ASL") = "m";

flight = timerange(eventfilter("LAUNCH"), eventfilter("APOGEE"), "closed");
openrocket = openrocket(flight, :);
burnout = openrocket(eventfilter("BURNOUT"), :);

% Simulation settings
% Assign every parameter to a single structure to use with struct2input() later
params.t_0 = seconds(burnout.Time);
params.dt = 0.01;
params.t_f = 40;

params.x_east_0 = burnout.("Lateral distance"); 
params.x_up_0 = burnout.("ASL");
params.v_east_0 = burnout.("Lateral velocity");
params.v_up_0 = burnout.("Vertical velocity");

% Vehicle settings
params.g = burnout.("Gravitational acceleration");
params.A_ref = (1e-2)^2 * (burnout.("Reference area")); % [m^2] reference area (from cm^2)
params.M_dry = burnout.Mass; % [kg]

% C_D lookup table from OR flight
C_D_lookup = openrocket(:, ["Mach number", "Drag coefficient"]);
C_D_lookup = sortrows(C_D_lookup, "Mach number", "ascend"); % Sort ascending by Mach number 
params.C_D_lookup = C_D_lookup;

params.C_D_brake = 1.2; % [-] approximately a flat plate
params.W_brake = 7.62 / 1e2; % [m] airbrake width - 3in (half dia.) 
params.N_brake = 2; % [-] how many airbrake fins there are

%% Create simulation cases
% Airbrake extension lengths
cases = table;
cases.L_brake = [0; 0.1; 0.2; 0.5; 1; 2; 3; 4; 5] / 1e2; % [cm]
disp(cases);

%% Run simulation cases
slx_name = "sim_2dof";
simin = struct2input(slx_name, params); % Put constants into simulation input
simins = table2inputs(simin, cases); % Create input array
simouts = sim(simins, ShowProgress = "on"); % simulate using input array

% Process each output
for si = 1:length(simouts)
    % store entire output in table cell
    % note {si} to create *cell array* to store 1 table in each row
    cases.output{si} = extractTimetable(simouts(si).yout);

    % summary stats
    cases.apogee(si) = max(cases.output{si}.x_up);
end
disp(cases);


%% Plots
figure(name = "Airbrake performance");
tiledlayout(2,1);

nexttile;
for ci = 1:height(cases)
    output = cases.output{ci};
    line_name = sprintf("L = %.f mm", 1e3 * cases.L_brake(ci));
    plot(output.Time, output.x_up - SL, DisplayName = line_name); % convert ASL to AGL
end
yline(3048, "--k", "10k Target", HandleVisibility = "off"); % so it doesn't show up on legend
ylabel("Altitude [m AGL]");
xlabel("Time");
legend(Location = "best");

nexttile;
plot(cases.L_brake * 1e3, cases.apogee - SL, "^k"); % convert M to mm, ASL to AGL
yline(3048, "--k", "10k Target", HandleVisibility = "off"); % so it doesn't show up on legend
ylabel("Apogee [m AGL]");
xlabel("Airbrake extension [mm]");
