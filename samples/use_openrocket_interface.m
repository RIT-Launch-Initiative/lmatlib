clear; close all;

%% Basic plots

otis = openrocket("data/OTIS.ork");
sim = otis.sims(1); % get simulation by number (name also works)
otis.simulate(sim); % execute simulation
data = openrocket.get_data(sim); % get all of the simulation's outputs

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
sim = otis.sims("20MPH-SA");
fins = otis.component("Trapezoidal Fins");
fh = fins.getHeight();

heights = fh * (0.8:0.05:1.2); % vary height from 80% to 120%

% Initialize output variables
ssm_launchrod = NaN(size(heights));
ssm_reference = NaN(size(heights));
ssm_burnout = NaN(size(heights));

% Refernce flight condition
ref_mach = 0.3;
ref_aoa = deg2rad(5);
ref_fcond = otis.flight_condition(ref_mach, ref_aoa);

for i = 1:length(heights)
    fins.setHeight(heights(i));
    
    ssm_reference(i) = otis.stability("LAUNCH", ref_fcond);

    % Simulate 
    data = otis.simulate(sim, stop = "apogee", outputs = "Stability margin"); 

    % Cut data to range of interest
    data_range = timerange(eventfilter("LAUNCHROD"), eventfilter("BURNOUT"), "openleft");
    data = data(data_range, :);
    % We need to data point immediately after the LAUNCHROD event because
    % OpenRocket only starts calculating stability margin after, not at, that
    % event. The "openleft" option for TIMERANGE cuts the table to include the
    % point immediately after the LAUNCHROD event. 

    ssm_launchrod(i) = data{1, "Stability margin"};
    ssm_burnout(i) = data{end, "Stability margin"};
end

figure(name = "Fin height sweep");
plot(100*heights, ssm_launchrod, DisplayName = "Rod clearance");
plot(100*heights, ssm_burnout, DisplayName = "Burnout");
plot(100*heights, ssm_reference, ...
    DisplayName = sprintf("Reference at Ma=%.1f, \\alpha=%.1f deg", ref_mach, rad2deg(ref_aoa)));
xlabel("Fin height [cm]"); ylabel("Stability [cal]");
title("Stability for range of fin heights, " + string(sim.getName))
legend;

yregion([-Inf 1.5], FaceColor = [0.8 0.1 0.1], HandleVisibility = "off");
yregion([4 Inf], FaceColor = [0.8 0.1 0.1], HandleVisibility = "off");

%% Monte Carlo simulations
n_sims = 20;
launch_spread = 4;
wind_speed_spread = 5;
wind_angle_spread = 60;

otis = openrocket("data/OTIS.ork");
sim = otis.sims("15MPH-SA");
opts = sim.getOptions();
opts.setWindTurbulenceIntensity(0.15);
opts.setLaunchRodDirection(deg2rad(60));

launch_angle = opts.getLaunchRodAngle();
wind_speed = opts.getWindSpeedAverage();

data_out = cell(1, n_sims);
for i = 1:length(data_out)
    opts.setLaunchRodAngle(deg2rad(launch_angle + (rand()-0.5)*launch_spread));
    opts.setLaunchIntoWind(false);
    opts.setWindDirection(-opts.getLaunchRodDirection + ...
        deg2rad((rand()-0.5)*wind_angle_spread));
    opts.setWindSpeedAverage(wind_speed + (rand()-0.5)*wind_speed_spread);

    data_out{i} = otis.simulate(sim, outputs = "ALL");
end

figure(name = "Landing point spread");
plot(0, 0, "+k", LineWidth = 2);
for i = 1:length(data_out)
    data = data_out{i};
    hit_data = data(eventfilter("GROUND_HIT"), :);
    plot(hit_data.("Position East of launch"), hit_data.("Position North of launch"), "xk");
end
xlabel("East [m]"); ylabel("North [m]");
title(sprintf("Landing point spread for %s\n Rod \\pm %.1f deg, Wind \\pm %.1f m/s at \\pm %.0f deg", ...
    sim.getName(), launch_spread/2, wind_speed_spread/2, wind_angle_spread/2))

%% Adjustable weight optimization
% Don't actually do this - we don't know any conditions to enough precision to
% make this number useful. But, it serves as an example.

otis = openrocket("data/OTIS.ork");
sim = otis.sims("15MPH-SA");
weight_opt_cost = make_cost_function(otis, sim, "Adjustable stability weight");

% Tolerance of 0.1m in altitude - more precise is not reasonable
% Tolerance of 10g in mass - more precise is not reasonable
% For fminsearch, tolerances are absolute - see https://www.mathworks.com/help/optim/ug/tolerance-details.html
opts = optimset(Display = "iter", TolFun = 0.1, TolX = 10e-3); 
opt_weight = fminsearch(weight_opt_cost, 1.5, opts);

% sim.getOptions().setWindTurbulenceIntensity(0);
otis.component("Adjustable stability weight").setComponentMass(opt_weight);
data = otis.simulate("15MPH-SA", outputs = "ALL");

figure(name = "Optimized trajectory");
plot_openrocket(data, "Altitude", "Vertical velocity", labels = ["BURNOUT", "APOGEE", "MAIN"]);
title(sprintf("Stability mass optimized to %.2f kg", opt_weight));


function func = make_cost_function(doc, sim, tuned_name)
    % All of these variables need to be pre-populated in the cost function, but
    % can't be arguments They can be made global, or the cost function is
    % nested within a function in which the variables are evaluated once.
    % 
    % https://www.mathworks.com/help/matlab/math/parameterizing-functions.html

    sim.getOptions().setWindTurbulenceIntensity(0);

    tuned_mass = doc.component(tuned_name);
    apogee_target = 3048;
    
    func = @cost;
    function f = cost(x)
        mass = x(1);
        tuned_mass.setComponentMass(mass);

        data = doc.simulate(sim, stop = "apogee", output = "Altitude");
        apogee_error = data{eventfilter("APOGEE"), "Altitude"} - apogee_target;
        f = abs(apogee_error);
    end        
end

