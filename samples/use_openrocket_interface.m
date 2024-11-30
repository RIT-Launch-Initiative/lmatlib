clear; close all;

%% Basic plots
otis = openrocket("data/OTIS.ork");
sim = otis.sims(1);
data = openrocket.get_data(sim);

figure(name = "Basic plots");
tiledlayout(2,1);

nexttile;
plot_openrocket(data, "Altitude", "Total velocity", ...
    end_ev = "GROUND_HIT", labels = ["BURNOUT", "APOGEE", "MAIN"]);

nexttile;
plot_openrocket(data, "Stability margin", "Angle of attack", ...
    start_ev = "LAUNCHROD", end_ev = "APOGEE", labels = ["LAUNCHROD", "BURNOUT"]);

%% Fin height sweep
otis = openrocket("data/OTIS.ork");
sim = otis.sims(1);
fins = otis.component("Trapezoidal Fins");
fh = fins.getHeight();

heights = fh * (0.8:0.05:1.2); % vary height from 80 to 120%
ssms = zeros(1, length(heights));
for i = 1:length(heights)
    fins.setHeight(heights(i));
    data = openrocket.get_data(sim, "Stability margin");
    ssms(i) = data{eventfilter("BURNOUT"), "Stability margin"};
end

figure(name = "Fin height sweep");
plot(100*heights, ssms);
xlabel("Fin height [cm]");
ylabel("Stability at BURNOUT [cal]");
yregion([4 Inf], FaceColor = [0.8 0.1 0.1]);

%% Wind speed sweep
otis = openrocket("data/OTIS.ork");
sim = otis.sims(1);
opts = sim.getOptions();

winds = 1:10;
aoa = zeros(1, length(winds));
for i = 1:length(winds)
    opts.setWindSpeedAverage(winds(i));
    data = openrocket.get_data(sim, "Angle of attack");
    aoa(i) = data{eventfilter("LAUNCHROD"), "Angle of attack"};
end

figure(name = "Wind speed sweep");
plot(winds, rad2deg(aoa));
xlabel("Average wind speed [m/s]");
ylabel("Angle of attack at LAUNCHROD [deg]");
yregion([15 Inf], FaceColor = [0.8 0.1 0.1]);

profile viewer;
