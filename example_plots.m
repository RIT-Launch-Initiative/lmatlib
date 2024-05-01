clear; close all;
import lplot.*;
load("sample_data/grim_test.mat");
load("sample_data/omen.mat");

%% Grim plots
figure(name = "pressure");
timeplot(grim_slow_t, "Pressure");
timeplot(rrc3_t, "Pressure");
legend("Grim", "RRC3");

figure(name = "temperature");
timeplot(grim_slow_t, "Temperature");
timeplot(rrc3_t, "Temperature");
legend("Grim", "RRC3");

descent = timerange(seconds(24.5), seconds(120));
figure(name = "altitude");
% a different way (I think better) way to name lines
timeplot(grim_slow_t, "Altitude", DisplayName = "Grim");
timeplot(rrc3_t, "Altitude", DisplayName = "RRC3");
timeplot(rs41_t(descent, :), "Altitude", "^r", DisplayName = "RS41"); % just the descent time interval
legend; % we defined them in DisplayName, so this just shows the legend

%% Calculate summary statistics and errors
fprintf("Apogee measured by RRC3: %.1f [m AGL]\n", max(rrc3_t.Altitude));
fprintf("Apogee measured by Grim: %.1f [m AGL]\n", max(grim_slow_t.Altitude));

% The RS-41 has the fewest data points, so interpolate to it in order to insert the least fake data
rs41_t_descent = rs41_t(descent, :);
rrc3_t_descent = retime(rrc3_t, rs41_t_descent.Time, "linear");
grim_t_descent = retime(grim_slow_t, rs41_t_descent.Time, "linear");

fprintf("RMS between Grim and RRC3: %.1f m\n", rms(grim_t_descent.Altitude - rrc3_t_descent.Altitude));
fprintf("RMS between RS-41 and RRC3: %.1f m\n", rms(rs41_t_descent.Altitude - rrc3_t_descent.Altitude));

%% OpenRocket flight report plots
% Multiple stacked time plots
figure(name = "OMEN flight report");
lay = tiledlayout(2,1);
nexttile;
plot_interval(omen, "Total velocity", "Vertical acceleration", ...
start_ev = "LAUNCHROD", end_ev = "APOGEE", labels = "BURNOUT");

nexttile;
plot_interval(omen, "Stability margin", "Angle of attack", ...
start_ev = "LAUNCHROD", end_ev = "APOGEE", labels = "BURNOUT");

stack_axes(lay);
