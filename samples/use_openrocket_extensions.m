%% Baseline
omen = openrocket("data/OMEN.ork");
sim = omen.sims("MATLAB");
site = launchsites("spaceport-america");

baseline_flight_data = omen.simulate(sim, outputs = "ALL");

baseline_drag_data = table;
baseline_drag_data.MACH = (0:0.1:2)';

for i_mach = 1:height(baseline_drag_data)
    fc = omen.flight_condition(baseline_drag_data.MACH(i_mach));
    [~, baseline_drag_data.DRAG(i_mach), ~, ~, ~] = ...
        omen.aerodata3(fc);
end

[imag, rast] = flight_basemap(site.lat, site.lon, 1e3);

figure(name = "Trajectory comparisons");
traj_ax = axes;
hold(traj_ax, "on");
grid(traj_ax, "on");
mapshow(imag, rast);
plot3_openrocket(traj_ax, baseline_flight_data, LineWidth = 2, DisplayName = "Baseline");
legend;

% List available models
disp("Available NCEP models:");
ncep.list
disp("Available GFS output grids: ");
ncep.list("gfs")
% Atmospheric model
times = datetime(2024, 06, 21, TimeZone = "MST") + hours([10 12]);
launchtime = datetime(2024, 06, 21, 10, 21, 00, TimeZone = "MST");

%% Custom atmosphere & wind models
% Get air data table

airdata = atmosphere("gfs", "pgrb2.1p00", site.lat, site.lon, launchtime, ...
    minpres = 450); % 400 mbar gives up 
% Celcius to Kelvin
airdata.TMP = airdata.TMP + 273.15;

atmos_flight_data = omen.simulate(sim, outputs = "ALL", ...
    atmos = airdata(:, ["HGT", "PRES", "TMP"]));
plot3_openrocket(traj_ax, atmos_flight_data, DisplayName = "Custom atmosphere");

wind_flight_data = omen.simulate(sim, outputs = "ALL", ...
    wind = airdata(:, ["HGT", "UGRD", "VGRD"]));
plot3_openrocket(traj_ax, wind_flight_data, DisplayName = "Custom wind");

models = import_rasaero_aerodata("data/OMEN_RA_Aerodata.csv");
drag_data = table;
drag_data.MACH = models.mach;
drag_data.DRAG = models.pick{"field", "CD", "aoa", 0}; % select raw values from xarray

drag_flight_data = omen.simulate(sim, outputs = "ALL", drag = drag_data);

plot3_openrocket(traj_ax, drag_flight_data, DisplayName = "Modified drag");

figure(name = "Drag lookup comparison");
hold on; grid on;
plot(baseline_drag_data.MACH, baseline_drag_data.DRAG, ...
    LineWidth = 2, DisplayName = "OpenRocket");
plot(drag_data.MACH, drag_data.DRAG, DisplayName = "RasAero II");
xlabel("Mach number");
ylabel("Drag coefficient");
legend;
xlim([0 2]);

figure(name = "Drag history")
hold on; grid on;
plot(baseline_flight_data.Time, baseline_flight_data.("Drag force"), ...
    LineWidth = 2, DisplayName = "OpenRocket");
plot(drag_flight_data.Time, drag_flight_data.("Drag force"), ...
    DisplayName = "RasAero II");
xlabel("Time");
ylabel("Drag coefficient");
legend;


