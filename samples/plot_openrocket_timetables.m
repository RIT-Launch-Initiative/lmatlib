% Examples for importing OpenRocket simulation data and plotting using built-in
% functions and the library's utilities
set(groot, "DefaultAxesNextPlot", "add");
clear; close all;

OR_source = "data/OMEN_OR_Output.csv";
simdata = import_openrocket_csv(OR_source);
set_grids("on");
set_interpreters("remove");

% Plot using native function
figure(name = "Mach number");
plot(simdata.Time, simdata.("Mach number"), DisplayName = "Mach number");
yregion([0.8 1.2], DisplayName = "Transonic region")
legend; % DisplayName arguments name each line as it is created
ylabel("Mach number [-]");
xlabel("Time");

% use event filter to pick out specific points in flight
apogee = eventfilter("APOGEE");
apogee_data = simdata(apogee, ["Altitude", "Total velocity"]);
fprintf("Apogee altitude (AGL): %.f m\n", apogee_data.Altitude);
fprintf("Velocity at apogee: %.1f m/s\n", apogee_data.("Total velocity"));

% Plot using utility
figure(name = "Stability");
tiledlayout(2,1)

nexttile;
plot_openrocket(simdata, "Altitude", ...
    labels = ["BURNOUT", "APOGEE", "MAIN"]);

nexttile; 
hold on;
plot_openrocket(simdata, "Stability margin", ...
    start_ev = "LAUNCHROD", end_ev = "APOGEE", labels = "BURNOUT");
yregion([-Inf 1], FaceColor = "#A2142F"); % too-low stability
yregion([4 Inf], FaceColor = "#A2142F"); % too-high stability
