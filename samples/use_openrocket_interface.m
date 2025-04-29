clear; close all;

run_sweep = true;
run_monte = true;
run_opt = true;

omen_path = "data/OMEN.ork"; %pfullfile("samples", "data", "OMEN.ork");

%% Basic plots
omen = openrocket(omen_path);
sim = omen.sims(1); % get simulation by number

drogue = omen.component(name = "Streamer"); % get streamer 
event = openrocket.get_deploy(drogue, sim); % get event for drogue chute
event.setDeployDelay(3); % 3-second drogue delay

openrocket.simulate(sim); % execute simulation
data = openrocket.get_data(sim); % get all of the simulation's outputs
% equivalently, data = openrocket.simulate(sim, outputs = "ALL") will do the same thing

figure(name = "Basic plots");
tiledlayout(2,1);

nexttile;
plot_openrocket(data, "Altitude", "Total velocity", ...
    end_ev = "GROUND_HIT", labels = ["BURNOUT", "APOGEE", "MAIN"]);

nexttile;
plot_openrocket(data, "Stability margin", "Angle of attack", ...
    start_ev = "LAUNCHROD", end_ev = "APOGEE", labels = ["LAUNCHROD", "BURNOUT"]);

drawnow; % so you get plots before waiting for the rest of this script to run
%% Fin height sweep
if run_sweep
    omen = openrocket(omen_path);
    sim = omen.sims("MATLAB");
    fins = omen.component(class = "FinSet"); 
    if ~isscalar(fins)
        error("Multiple fin sets found");
    end
    fh = fins.getHeight();

    heights = fh * (0.8:0.05:1.2); % vary height from 80% to 120%

    % Initialize output variables
    ssm_launchrod = NaN(size(heights));
    ssm_reference = NaN(size(heights));
    ssm_burnout = NaN(size(heights));

    % Refernce flight condition
    ref_mach = 0.3;
    ref_aoa = deg2rad(5);
    ref_fcond = omen.flight_condition(ref_mach, ref_aoa);

    parfor i_sim = 1:length(heights);
    
        pardoc = openrocket(omen_path);
        parsimobj = pardoc.sims("MATLAB")
        fins = pardoc.component(class = "FinSet"); 

        fins.setHeight(heights(i_sim));
        
        % ssm_reference(i_sim) = pardoc.stability("LAUNCH", ref_fcond);

        % Simulate 
        data = openrocket.simulate(parsimobj, outputs = "Stability margin"); 

        % Cut data to range of interest
        data_range = timerange(eventfilter("LAUNCHROD"), eventfilter("BURNOUT"), "openleft");
        data = data(data_range, :);
        % We need to data point immediately after the LAUNCHROD event because
        % OpenRocket only starts calculating stability margin after, not at, that
        % event. The "openleft" option for TIMERANGE cuts the table to include the
        % point immediately after the LAUNCHROD event. 

        ssm_launchrod(i_sim) = data{1, "Stability margin"};
        ssm_burnout(i_sim) = data{end, "Stability margin"};
    end

    figure(name = "Fin height sweep");
    hold on; grid on;
    plot(100*heights, ssm_launchrod, DisplayName = "Rod clearance");
    plot(100*heights, ssm_burnout, DisplayName = "Burnout");
    plot(100*heights, ssm_reference, ...
        DisplayName = sprintf("Reference at Ma=%.1f, \\alpha=%.1f deg", ref_mach, rad2deg(ref_aoa)));
    xlabel("Fin height [cm]"); ylabel("Stability [cal]");
    title("Stability for range of fin heights, " + string(sim.getName))
    legend;

    yregion([-Inf 1.5], FaceColor = [0.8 0.1 0.1], HandleVisibility = "off");
    yregion([4 Inf], FaceColor = [0.8 0.1 0.1], HandleVisibility = "off");
    drawnow; % so you get plots before waiting for the rest of this script to run
end

%% Monte Carlo simulations
if run_monte
    n_sims = 20;
    launch_spread = 4;
    wind_speed_spread = 5;
    wind_angle_spread = 60;

    omen = openrocket(omen_path);
    sim = omen.sims("MATLAB");
    opts = sim.getOptions();
    opts.setWindTurbulenceIntensity(0.15);
    opts.setLaunchRodDirection(deg2rad(60));
    opts.setTimeStep(0.05); % lower rate to improve performacne

    launch_angle = opts.getLaunchRodAngle();
    wind_speed = opts.getWindSpeedAverage();

    data_out = cell(1, n_sims);
    for i = 1:length(data_out)
        opts.setLaunchRodAngle(deg2rad(launch_angle + (rand()-0.5)*launch_spread));
        opts.setLaunchIntoWind(false);
        opts.setWindDirection(-opts.getLaunchRodDirection + ...
            deg2rad((rand()-0.5)*wind_angle_spread));
        opts.setWindSpeedAverage(wind_speed + (rand()-0.5)*wind_speed_spread);

        data_out{i} = openrocket.simulate(sim, outputs = "ALL");
    end

    figure(name = "Landing point spread");
    hold on;
    plot(0, 0, "+k", LineWidth = 2);
    for i = 1:length(data_out)
        data = data_out{i};
        hit_data = data(eventfilter("GROUND_HIT"), :);
        plot(hit_data.("Position East of launch"), hit_data.("Position North of launch"), "xk");
    end
    xlabel("East [m]"); ylabel("North [m]");
    title(sprintf("Landing point spread for %s\n Rod \\pm %.1f deg, Wind \\pm %.1f m/s at \\pm %.0f deg", ...
        sim.getName(), launch_spread/2, wind_speed_spread/2, wind_angle_spread/2))
end

%% Adjustable weight optimization
if run_opt
    omen = openrocket(omen_path);
    sim = omen.sims("MATLAB");
    opt_mass_name = "Mass Component";
    apogee_target = 3048; % [m] 
    apogee_tolerance = 1; % [m] stop sim when we are within this value of the desired apogee

    % see https://www.mathworks.com/help/optim/ug/output-function.html#f11454
    weight_opt_cost = make_cost_function(omen, sim, opt_mass_name, 3048);
    stop_close_enough = @(~, optim_values, ~) optim_values.fval <= apogee_tolerance; 
    opts = optimset(Display = "iter", OutputFcn = stop_close_enough); 
    opt_weight = fminsearch(weight_opt_cost, 1.5, opts);

    % sim.getOptions().setWindTurbulenceIntensity(0);
    omen.component(name = opt_mass_name).setComponentMass(opt_weight);
    data = openrocket.simulate(sim, outputs = "ALL");

    figure(name = "Optimized trajectory");
    plot_openrocket(data, "Altitude", "Vertical velocity", labels = ["BURNOUT", "APOGEE", "MAIN"]);
    title(sprintf("Stability mass optimized to %.2f kg", opt_weight));
    drawnow; % so you get plots before waiting for the rest of this script to run
end

function func = make_cost_function(doc, sim, tuned_name, apogee_target)
    % All of these variables need to be pre-populated in the cost function, but
    % can't be arguments They can be made global, or the cost function is
    % nested within a function in which the variables are evaluated once.
    % 
    % https://www.mathworks.com/help/matlab/math/parameterizing-functions.html

    sim.getOptions().setWindTurbulenceIntensity(0);
    tuned_mass = doc.component(name = tuned_name);
    func = @cost;
    function f = cost(x)
        mass = x(1);
        tuned_mass.setComponentMass(mass);

        data = doc.simulate(sim, outputs = "Altitude");
        apogee_error = data{eventfilter("APOGEE"), "Altitude"} - apogee_target;
        f = abs(apogee_error);
    end        
end
