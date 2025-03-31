clear; close all;

run_sweep = false;
run_monte = false;
run_opt = false;
run_atmos = true;
run_wind = true;
run_drag = true;

otis_path = "data/OTIS.ork"; %pfullfile("samples", "data", "OTIS.ork");

%% Basic plots
otis = openrocket(otis_path);
sim = otis.sims(1); % get simulation by number

drogue = otis.component(name = "Streamer"); % get streamer 
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


%% Fin height sweep
if run_sweep
    otis = openrocket(otis_path);
    sim = otis.sims("20MPH-SA");
    fins = otis.component(class = "FinSet"); 
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
    ref_fcond = otis.flight_condition(ref_mach, ref_aoa);

    for i = 1:length(heights)
        fins.setHeight(heights(i));
        
        ssm_reference(i) = otis.stability("LAUNCH", ref_fcond);

        % Simulate 
        data = openrocket.simulate(sim, stop = seconds(3), outputs = "Stability margin"); 
        % stop at 3 seconds to not simulate much after burnout - significant performance benefit

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
end

%% Monte Carlo simulations
if run_opt
    n_sims = 20;
    launch_spread = 4;
    wind_speed_spread = 5;
    wind_angle_spread = 60;

    otis = openrocket(otis_path);
    sim = otis.sims("15MPH-SA");
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
    % Don't actually do this - we don't know any conditions to enough precision to
    % make this number useful. But, it serves as an example.

    otis = openrocket(otis_path);
    sim = otis.sims("15MPH-SA");
    weight_opt_cost = make_cost_function(otis, sim, "Adjustable stability weight");

    % see https://www.mathworks.com/help/optim/ug/output-function.html#f11454
    stop_close_enough = @(~, optim_values, ~) optim_values.fval <= 0.1; % stop sim when we are within 0.1m
    opts = optimset(Display = "iter", OutputFcn = stop_close_enough); 
    opt_weight = fminsearch(weight_opt_cost, 1.5, opts);

    % sim.getOptions().setWindTurbulenceIntensity(0);
    otis.component(name = "Adjustable stability weight").setComponentMass(opt_weight);
    data = openrocket.simulate(sim, outputs = "ALL");

    figure(name = "Optimized trajectory");
    plot_openrocket(data, "Altitude", "Vertical velocity", labels = ["BURNOUT", "APOGEE", "MAIN"]);
    title(sprintf("Stability mass optimized to %.2f kg", opt_weight));
end

%% Custom Simulation Inputs

% Baseline
if run_atmos || run_wind || run_drag
    otis = openrocket(otis_path);
    sim = otis.sims("20MPH-SA");

    baseline_flight_data = otis.simulate(sim, outputs = "ALL");

    baseline_drag_data = table;
    baseline_drag_data.MACH = (0:0.1:2)';

    for i_mach = 1:height(baseline_drag_data)
        fc = otis.flight_condition(baseline_drag_data.MACH(i_mach));
        [~, baseline_drag_data.DRAG(i_mach), ~, ~, ~] = ...
            otis.aerodata3(fc);
    end
    
    figure(name = "Trajectory comparisons");
    traj_ax = axes;
    hold(traj_ax, "on");
    grid(traj_ax, "on");
    plot_trajectory(traj_ax, baseline_flight_data, LineWidth = 2, DisplayName = "Baseline");
    legend;
end

% Atmospheric model
if run_atmos
    site = launchsites("spaceport-america");
    times = datetime(2024, 06, 21, TimeZone = "MST") + hours([10 12]);
    launchtime = datetime(2024, 06, 21, 10, 21, 00, TimeZone = "MST");
    refs = ncep.anl("gfs", "pgrb2.1p00", times);
    refs.attach("data");
    atmos = atmosphere.from_ncep(refs, lats = site.lat + [-1 1], ...
        lons = site.lon + [-1 1]);
    ac = atmos.aircolumn(site.lat, site.lon, launchtime);
    airdata = table;
    airdata.HGT = ac.pick{"field", "HGT"};
    airdata.PRES = 100*str2double(extract(ac.layer, digitsPattern));
    airdata.TMP = ac.pick{"field", "TMP"} + 273.15;

    atmos_flight_data = otis.simulate(sim, outputs = "ALL", atmos = airdata);
    plot_trajectory(traj_ax, atmos_flight_data, DisplayName = "Custom atmosphere");
end

% Wind model
if run_wind
    site = launchsites("spaceport-america");
    times = datetime(2024, 06, 21, TimeZone = "MST") + hours([10 12]);
    launchtime = datetime(2024, 06, 21, 10, 21, 00, TimeZone = "MST");
    refs = ncep.anl("gfs", "pgrb2.1p00", times);
    refs.attach("data");
    atmos = atmosphere.from_ncep(refs, lats = site.lat + [-1 1], ...
        lons = site.lon + [-1 1]);
    ac = atmos.aircolumn(site.lat, site.lon, launchtime);
    winddata = table;
    winddata.HGT = ac.pick{"field", "HGT"};
    winddata.UGRD = ac.pick{"field", "UGRD"};
    winddata.VGRD = ac.pick{"field", "VGRD"};

    wind_flight_data = otis.simulate(sim, outputs = "ALL", wind = winddata);
    plot_trajectory(traj_ax, wind_flight_data, DisplayName = "Custom wind");
end

% Drag model
if run_drag
    models = import_rasaero_aerodata("data/OMEN_RA_Aerodata.csv");
    drag_data = table;
    drag_data.MACH = models.mach;
    drag_data.DRAG = models.pick{"field", "CD", "aoa", 0};
    
    drag_flight_data = otis.simulate(sim, outputs = "ALL", drag = drag_data);

    plot_trajectory(traj_ax, drag_flight_data, DisplayName = "Modified drag");

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

end


function func = make_cost_function(doc, sim, tuned_name)
    % All of these variables need to be pre-populated in the cost function, but
    % can't be arguments They can be made global, or the cost function is
    % nested within a function in which the variables are evaluated once.
    % 
    % https://www.mathworks.com/help/matlab/math/parameterizing-functions.html

    sim.getOptions().setWindTurbulenceIntensity(0);

    tuned_mass = doc.component(name = tuned_name);
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

function ph = plot_trajectory(ax, flight_data, varargin)
    ph = plot3(ax, flight_data.("Position East of launch"), ....
        flight_data.("Position North of launch"), ...
        flight_data.Altitude, varargin{:});

    xlabel(ax, "East [m]");
    ylabel(ax, "North [m]");
    zlabel(ax, "Altitude [m]");
    % daspect([1 1 1]);
end
