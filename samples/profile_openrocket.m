clear; close all;

profile on;
%% Simulate entire time
otis = openrocket("data/OTIS.ork");
sim = otis.sims(1);
opts = sim.getOptions();
opts.setTimeStep(0.05);

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

%% Up to apogee
otis = openrocket("data/OTIS.ork");
sim = otis.sims(1);
opts = sim.getOptions();
opts.setTimeStep(0.05);

winds = 1:10;
aoa = zeros(1, length(winds));
for i = 1:length(winds)
    opts.setWindSpeedAverage(winds(i));
    openrocket.simulate(sim, "bypass", "apogee");
    data = openrocket.get_data(sim, "Angle of attack");
    aoa(i) = data{eventfilter("LAUNCHROD"), "Angle of attack"};
end

figure(name = "Wind speed sweep");
plot(winds, rad2deg(aoa));
xlabel("Average wind speed [m/s]");
ylabel("Angle of attack at LAUNCHROD [deg]");
yregion([15 Inf], FaceColor = [0.8 0.1 0.1]);
profile viewer;
