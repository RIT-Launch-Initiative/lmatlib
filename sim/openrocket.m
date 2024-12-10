%% OPENROCKET
% MATLAB utility class that loads an OpenRocket application and wraps
% underlying Java methods for automating OpenRocket calculations and
% simulations.
% 
% OVERVIEW
% The class provides the user with the ability to load an OpenRocket program
% instance, open a document, and save it. Minimally, the `document` field is
% sufficient to do everything else. All of the objects are reference-type
% (handle), not copy-on-write as you may be used to in MATLAB. Unless you
% deep-copy an object, modifying members of this class will modify the
% underlying OpenRocket document. Very common operations, namely running
% simulations, finding components, and basic inertial and aerodynamic outputs
% are also wrapped into helper functions for ease-of-use. 
% 
% If you are programmatically modifying the simulation, you will likely need to
% use the relevant Java object's methods directly (getOptions(), setHeight(),
% ...). They are usually self-explanatory (if somewhat cumbersome).
% Tab-completion works for Java object methods, but to get a complete list of
% the methods available, you can use methodsview(class(java_object)), for example:
%   or = openrocket("test.ork");
%   fins = or.shortname("fins");
%   methodsview(class(fins))
%
% In some cases, parts of the rocket or simulation are held under
% FlightConfigurations. Quantities dependent on FlightConfigurations, such as
% aerodynamics, ineritas, or recovery events, will require a
% FlightConfiguration to operate on. Usually, if one is not provided, they
% assume the selected configuration is the desired one.
%   
% PREREQUISITES
%   - OpenRocket installation, containing jre/ and OpenRocket.jar class file
%   - Version of MATLAB <a href="https://www.mathworks.com/support/requirements/openjdk.html">compatible</a> with the Java version used by OpenRocket. 
%   At time of writing, the latest release is OR 23.09 using JDK 17, which
%   requires MATLAB 2024a or newer
%   - Set up Java class path and Java environment using openrocket_setup(...)
% 
%% LISTING
% Static utilities
% openrocket.get_data(sim[, variables])         :retrieve data from up-to-date simulation
% openrocket.search_components(rkt, mode, key)  :search Rocket object's
%   children for components by name or class
% openrocket.list_variables()                   :list possible flight data variables
% -----------------------------------------------------------------------------------
% 
% Document methods
% or = openrocket(path)                 :construct interface from path to .ork file
% or.save()                             :save modiied OpenRocket file to original document
% rkt = or.rocket()                     :return document's Rocket object
% cfg = or.get_config([identifier])     :get FlightConfiguration by name or number,
%   returning active one if unspecified
% or.set_config(cfg)                    :set active FlightConfiguration
% part = or.component(name)             :get component(s) by name
% part = or.shortname(name)             :get common key components by name -
%   "main", "drogue", etc.
% -----------------------------------------------------------------------------------
% Simulation methods
% 
% sim = or.sims([identifier])           :retrieve Simulation by name or number,
%   returning all if unspecified
% [data =] or.simulate(identifier[, stop = "apogee"][, outputs = [...]])
%   Execute simulation (using the same identifier as before), optionally
%   returning data through "outputs".
% -----------------------------------------------------------------------------------
% Aerodynamic and inertial calculations
% (Unless otherwise specified, these methods use the active FlightConfiguration)
% 
% fc = or.flight_condition([mach, aoa, theta, rpy_rate])
%   Get FlightCondition object to use for aerodynamic methods
% data = or.aerodata6(fc)       :get AerodynamicForces object for 6-DOF simulation
% [CP, CD, CN, Cm, CNa] = or.aerodata3(fc)
%   Get 3-DOF parameters
% [d, A] = obj.refdims()        :get reference length and area
% [CG, mass, moi] = obj.massdata(state)
%   Get inertial properties at LAUNCH or BURNOUT
% ssm = obj.stability(state, fc)
%   Get stability at LAUNCH or BURNOUT, at specified flight condition


% Written by Yevgeniy Gorbachev
% November 2024

%% NOTES FOR DEVELOPERS
% - The application startup is translated from <a href="https://github.com/SilentSys/orhelper/blob/master/orhelper/_orhelper.py">SilentSys/orhelper</a>, many thanks 
% for doing the annoying part.
% - There is no significant performance difference between <orhelper> and this
% class, because the bulk of the performance cost is the simulation itself.
% - Operations cannot be parallelized becuase PARFOR <a href="https://www.mathworks.com/help/parallel-computing/objects-and-handles-in-parfor-loops.html">requires</a> the data 
% passed in to support SAVE and LOAD, which Java objects in MATLAB do not
 
classdef (Sealed) openrocket < handle
    
    properties (SetAccess = private, GetAccess = public)
        document; % OpenRocketDocument Java object
    end

    properties (Constant, Access = private)
        started = openrocket.start(); % dummy constant to ensure start() is called once

        types = openrocket.make_type_map(); % Column name to FlightDataType
        units = openrocket.make_units_map(); % Column name to units

        saver = net.sf.openrocket.file.GeneralRocketSaver(); % Rocket file saver
        barrowman = net.sf.openrocket.aerodynamics.BarrowmanCalculator(); % Barrowman calculator
        masscalc = net.sf.openrocket.masscalc.MassCalculator();  % Mass calculator
        warnings = net.sf.openrocket.logging.WarningSet(); % Empty warning set (for Barrowman calculator input)
        status_uptodate = "UPTODATE";  % Magic value: Indicates up-to-date simulation
        time_name = "Time"; % Magic value: Name given to the time axis
        reco_name = "RECOVERY_DEVICE_DEPLOYMENT"; % Magic value: Name printed for recovery events
    end

    properties (Access = private)
        file;   % Java File object
        loader; % Java Loader object
    end

%% STATIC UTILITIES
    methods (Static, Access = public)

        function data = get_data(sim, variables)
            % Retreive data from an OpenRocket simulation
            % data = get_data(sim, variables)
            %   sim         Simulation object (returned by openrocket.sims() or directly accessed)
            %   variables   (Optional) list of variables to return
            %               Defaults to add all of them
            % 
            % NOTE Does not automatically simulate outdated simulation
            % NOTE All variables are unitless or MKS
            arguments
                sim
                variables string = openrocket.types.keys;
            end

            if string(sim.getStatus) ~= openrocket.status_uptodate
                error("Simulation not up to date. Please remember to run openrocket.simulate().");
            end

            variables = variables(variables ~= openrocket.time_name);
            branch = sim.getSimulatedData().getBranch(0); % NOTE assumes branch 0

            % Create destination timetable
            time = openrocket.jarr2double(branch.get(openrocket.types(openrocket.time_name)));
            data = timetable(Size = [length(time) length(variables)], ...
                VariableTypes = repmat("double", 1, length(variables)), ...
                VariableNames = variables, ...
                RowTimes = seconds(time)); % NOTE assumes time is in seconds

            % Assign data
            for idx = 1:length(variables)
                col = variables(idx);
                jarr = branch.get(openrocket.types(col));
                if isempty(jarr)
                    data = removevars(data, col);
                else
                    data.(col) = openrocket.jarr2double(branch.get(openrocket.types(col)));
                end
            end

            % Assign metadata
            evs = toArray(branch.getEvents());
            names = arrayfun(@(ev) string(ev.getType().name()), evs);
            times = arrayfun(@(ev) ev.getTime(), evs);
            reco_idx = find(names == openrocket.reco_name);
            switch length(reco_idx)
                case 1
                    names(reco_idx(1)) = "MAIN";
                case 2
                    names(reco_idx(1)) = "DROGUE";
                    names(reco_idx(2)) = "MAIN";
            end

            data.Properties.Events = eventtable(seconds(times), EventLabels = names); 
            data.Properties.VariableUnits = openrocket.units(data.Properties.VariableNames);
            data.Properties.VariableContinuity = repmat("continuous", 1, width(data));
            data.Properties.Description = "openrocket"; % Identify as OpenRocket import
        end

        function parts = search_components(rkt, mode, key)
            % Search rocket for components
            % parts = search_components(rkt, mode, key)
            %   rkt     Rocket object (or subcomponent)
            %   mode    Search by name or Java class
            %   key     If mode is "name", the component's (case-sensitive) editor name.
            %           If mode is "class", a (super) class under openrocket.rocketcomponent.(...)
            %   
            %   EXAMPLES
            %   or = openrocket("test.ork");
            %   fins = openrocket.search_components(or.rocket(), class = "FinSet");
            %   streamer = openrocket.search_components(or.rocket(), name = "Streamer");

            arguments
                rkt
                mode (1,1) string;
                key (1,1) string;
            end 

            if ~isa(rkt, "net.sf.openrocket.rocketcomponent.RocketComponent")
                error("First argument must be a Rocket (or component thereof)");
            end

            parts = rkt.getAllChildren().toArray();
            switch (mode)
                case "name"
                    match = string(parts) == key;
                case "class"
                    jclass = sprintf("net.sf.openrocket.rocketcomponent.%s", key);
                    match = false(size(parts));
                    for i_part = 1:numel(parts)
                        match(i_part) = isa(parts(i_part), jclass);
                    end
                otherwise
                    error("Search mode '%s' not recognized", mode);
            end
            parts = parts(match);
        end

        function vars = list_variables
            % List available flight data types
            vars = openrocket.types.keys;
        end
    end

%% COMPONENT METHODS
    methods 
        function obj = openrocket(orkpath)
            % Construct openrocket object from path to OpenRocket file
            arguments
                orkpath (1,1) string
            end

            [status, info] = fileattrib(orkpath);
            if ~status
                error("File '%s' does not exist", orkpath);
            end
            path = info.Name;
            [~, ~, ext] = fileparts(path);
            if ext ~= ".ork"
                error("File has extension '%s' instead of '.ork'", ext);
            end

            obj.file = java.io.File(path);
            obj.loader = net.sf.openrocket.file.GeneralRocketLoader(obj.file);
            obj.document = obj.loader.load();
        end

        function save(obj)
            % Save OpenRocket document to file
            obj.saver.save(obj.file, obj.document);
        end

        function out = sims(obj, sim)
            % Get simulations 
            %   sim = or.sims()       
            %       returns all simulations
            %   sim = or.sims(index)  
            %       return simulation at indices (1-based)
            %   sim = or.sims(name)   
            %       return all simulations matching name

            sims = toArray(obj.document.getSimulations());
            if nargin == 1
                out = sims; 
            elseif isnumeric(sim)
                out = sims(sim);
            elseif isstring(sim)
                names = arrayfun(@(sim) string(sim.getName()), sims);
                out = sims(names == sim);
            else 
                error("Identifier not supported");
            end
        end

        function data = simulate(obj, ident, params)
            % Simulate an OpenRocket simulation 
            %   or.simulate(ident, stop = "", outputs = "")
            %       ident       Simulation number, or name, or object (returned fom sims(...))
            %       stop        (Optional) Currently, only supports "apogee" to
            %                   stop sim at apogee
            %       outputs     (Optional) Convert flight data
            %                   String array e.g. ["Altitude", "Stability margin"]
            %                   OUT 
            %                   true (logical) outputs all values
            %   EXAMPLES
            %       or = openrocket("data/OTIS.ork");
            %       sim = or.sims(1);
            %       % Simulate using object, up to apogee
            %       or.sim(sim, stop = "apogee"); 
            %       % Simulate first simulation and assign all outputs
            %       data = or.sim(1, output = "ALL");  
            %       % Simulate named simulation, get stability data
            %       data = or.sim("15MPH-SA", output = "Stability margin"); 

            arguments
                obj openrocket
                ident
                params.stop (1, :) = [];
                params.outputs (1, :) string = [];
            end


            % NOTE if this is modified, make sure the input to simulate() is
            % exactly the class Java expects - if it isn't, MATLAB thinks
            % you're trying to call simulate() from some random Toolbox instead
            % of giving a Java error.

            import net.sf.openrocket.simulation.listeners.system.*
            import net.sf.openrocket.simulation.extension.example.*
            listener_class = "net.sf.openrocket.simulation.listeners.SimulationListener"; 
            
            if isa(ident, "net.sf.openrocket.document.Simulation")
                sim = ident;
            else
                sim = obj.sims(ident);
            end

            if length(sim) ~= 1
                error("More than one simulation specified:\n%s", mat2str(sims));
            end

            outs = params.outputs;
            if isempty(outs) && nargout > 0
                error("Simulation output assigned, but not requested through 'outputs = ...'");
            elseif ~isempty(outs) && nargout == 0
                warning("Simulation output requested through 'outputs = ....', but not assigned.");
            end


            listeners = javaArray(listener_class, 0); %#ok 
            % listeners has an #ok on both lines because the sim call is in an
            % evalc() so MATLAB doesn't know we're using this variable and warns
            extensions = java.util.List.of();
            if isstring(params.stop) && params.stop == "apogee"
                listeners = javaArray(listener_class, 1);
                listeners(1) = OptimumCoastListener.INSTANCE; %#ok
            elseif isduration(params.stop)
                stopper = StopSimulation();
                stopper.setStopStep(intmax);
                stopper.setReportRate(intmax);
                stopper.setStopTime(seconds(params.stop));
                extensions = java.util.List.of(stopper);
            end

            sim.copyExtensionsFrom(extensions);

            % The stop simulation listener feels the need to print directly to
            % the command line instead of the logging queue, so we eat sim's
            % output with evalc so it doesn't pollute our command line.
            evalc("sim.simulate(listeners)");

            sim.copyExtensionsFrom(java.util.List.of());
            
            if ~isempty(outs) && nargout > 0
                if outs == "ALL"
                    data = openrocket.get_data(sim);
                else
                    data = openrocket.get_data(sim, outs);
                end
            end
        end

        function rkt = rocket(obj)
            % Get Rocket object 
            %   rkt = or.rocket();
            rkt = obj.document.getRocket();
        end

        function cfg = get_config(obj, ident)
            % Get FlightConfiguration 
            %   config = or.get_config
            %       return selected flight configuration
            %   config = or.get_config(1)
            %       return indexed flight configuration
            %   config = or.get_config("[No motors]")
            %       return named flight configuration

            if nargin == 1
                cfg = obj.rocket().getSelectedConfiguration();
            elseif nargin == 2
                configs = toArray(obj.rocket().getFlightConfigurations());
                if isnumeric(ident)
                    cfg = configs(ident);
                elseif isstring(ident)
                    cfg = configs(string(configs) == ident);
                end
                if isempty(cfg)
                    error("No flight configuration found");
                end
            end
        end

        function set_config(obj, ident)
            % Set currently selected FlightConfiguration
            cfg = obj.get_config(ident);
            obj.rocket().setSelectedConfiguration(cfg.getId);
        end

        function parts = component(obj, name)
            % Search entire Rocket for component name
            %   component = or.component(name)

            rkt = obj.rocket();
            parts = openrocket.search_components(rkt, name = name);
            if isempty(parts)
                error("Could not find part '%s'", name);
            elseif length(parts) >= 2
                error("Found %d components named '%s'", length(parts), name);
            end
        end

        function part = shortname(obj, name, multiple)
            % Search Rocket for generic component
            %   part = or.shortname(name, multiple = false)
            %   name is one of
            %           "drogue": recovery device configured for APOGEE or EJECTION
            %           "main": recovery device configured for ALTITUDE, or the
            %           recovery device if there is only one, regardless of configuration.
            %           "fins": fin set
            %           "mount": motor mount
            %           "nose": nose cone
            %   multiple = true allows return of multiple parts, defaults to false
            %
            arguments
                obj openrocket;
                name (1,1) string;
                multiple (1,1) logical = false;
            end

            rkt = obj.rocket();
            cfg = obj.get_config();

            switch name
                % repeated code is ugly, but not repeating it will be uglier
                case "drogue" % all parachutes set to deploy at apogee or with ejection charge
                    chutes = openrocket.search_components(rkt, class = "RecoveryDevice");
                    to_return = false(size(chutes));
                    for i = 1:numel(chutes)
                        dep_cfg = chutes(i).getDeploymentConfigurations.get(cfg.getId);
                        ev = dep_cfg.getDeployEvent;
                        to_return(i) = ev == ev.APOGEE || ev == ev.EJECTION;
                        % APOGEE and EJECTION are members of <ev> becuase it is
                        % not straightforward to get at anonymous class members
                        % more directly
                    end
                    part = chutes(to_return);
                case "main" % all parachutes set to deploy at specific altitude
                    chutes = openrocket.search_components(rkt, class = "RecoveryDevice");
                    to_return = false(size(chutes));
                    for i = 1:numel(chutes)
                        dep_cfg = chutes(i).getDeploymentConfigurations.get(cfg.getId);
                        to_return(i) = dep_cfg.getDeployEvent == dep_cfg.getDeployEvent.ALTITUDE;
                    end
                    % If there is only one parachute, it is always main
                    if numel(chutes) == 1
                        to_return = true;
                    end
                    part = chutes(to_return);
                case "fins"
                    part = openrocket.search_components(rkt, class = "FinSet");
                case "mount"
                    mounts = openrocket.search_components(rkt, class = "MotorMount");
                    to_return = false(size(mounts));
                    for i = 1:numel(mounts)
                        to_return(i) = mounts(i).isMotorMount;
                    end
                    part = mounts(to_return);
                case "nose"
                    part = openrocket.search_components(rkt, class = "NoseCone");
                otherwise 
                    error("Short-name '%s' not recognized");
            end

            if isempty(part) 
                error("No part matching '%s' found");
            elseif ~multiple && ~isscalar(part)
                error("Multiple parts matching '%s' found", name);
            end
        end

    end


%% AERODYNAMIC CALCULATIONS
    methods
        function fc = flight_condition(obj, mach, aoa, theta, rpy_rate)
            % Get FlightCondition object for given inputs
            %   fc = flight_condition(obj, mach, aoa, theta, rpy_rate)
            %       returns object with properties populated in indicated sequence
            %       Some or all can remain unspecified, and default to 0 (except Ma=0.3)
            %   
            %   EXAMPLES
            %   fc = or.flight_condition() - default flight condition
            %   fc = or.flight_condition(0.2, deg2rad(5)) - Ma=0.2, alpha=5 deg

            arguments
                obj openrocket;
                mach (1,1) double = 0.3;
                aoa (1,1) double = 0;
                theta (1,1) double = 0;
                rpy_rate (3,1) double = [0; 0; 0];
            end
            cfg = obj.get_config();
            fc = net.sf.openrocket.aerodynamics.FlightConditions(cfg);
            fc.setMach(mach);
            fc.setAOA(aoa);
            fc.setTheta(theta);
            fc.setRollRate(rpy_rate(1));
            fc.setPitchRate(rpy_rate(2));
            fc.setYawRate(rpy_rate(3));
        end

        function [data] = aerodata6(obj, fc)
            % Get entire (6-DOF) AerodynamicForces object at given FlightCondition
            %   forces = or.aerodata6(fc)
            %       returns an AerodynamicForces object, in which individual
            %       coefficients are accessed through their Java accessors

            cfg = obj.get_config();
            data = openrocket.barrowman.getAerodynamicForces(cfg, fc, openrocket.warnings);
        end

        function [CP, CD, CN, Cm, CNa] = aerodata3(obj, fc)
            % Get 3-DOF aerodynamic data at given FlightCondition
            %   [CP, CD, CN, Cm, CNa] = or.aerodata3(fc)
            %       CP      Center of pressure (from nose, m)
            %       CD      Drag coefficient
            %       CN      Normal force coefficient
            %       CM      Pitching moment coefficient
            %       CNa     Normal force derivative

            cfg = obj.get_config();
            data = openrocket.barrowman.getAerodynamicForces(cfg, fc, openrocket.warnings);
            CP = data.getCP.x;
            CD = data.getCD;
            CN = data.getCN;
            Cm = data.getCm;
            CNa = data.getCNa;
        end

        function [d, A] = refdims(obj)
            % Get reference dimensions
            %   [d, A] = or.refdims()
            %       d   reference length (m)
            %       A   reference area (m^2)

            fc = obj.flight_condition();
            d = fc.getRefLength();
            A = fc.getRefArea();

        end

        function [CG, mass, moi] = massdata(obj, state)
            % Get mass data for rocket state in selected flight state
            %   [CG, mass, moi] = or.massdata(state)
            %       state   "LAUNCH" or "BURNOUT"
            %       
            %       CG      Center of mass (from nose, m)
            %       mass    Vehicle mass (kg)
            %       moi     Principal moments about [roll; pitch; yaw] (kg*m^2)

            cfg = obj.get_config();
            switch state
                case "LAUNCH"
                    data = openrocket.masscalc.calculateLaunch(cfg);
                case "BURNOUT"
                    data = openrocket.masscalc.calculateBurnout(cfg);
                otherwise 
                    error("Unrecognized flight state '%s'", state)
            end
            CG = data.getCM.x;
            mass = data.getMass;
            moi = [data.getIxx; data.getIyy; data.getIzz];
        end
        
        function ssm = stability(obj, state, fc)
            % Get stability margin given flight state and flight condition
            %   ssm = or.massdata(state, fc)
            %       state   "LAUNCH" or "BURNOUT"
            %       fc      FlightCondition
            %       
            %       ssm     Stability margin, in calibers

            cfg = obj.get_config();
            CP3 = openrocket.barrowman.getCP(cfg, fc, openrocket.warnings);
            CP = CP3.x;
            [CG, ~, ~] = obj.massdata(state);
            ssm = (CP - CG) / fc.getRefLength();
        end
    end
    
%% INTERNAL UTILITIES
    methods (Static, Access = private)
        function ret = start
            % Creates barebones OpenRocket application. 
            % This is intended to be called exactly once in a MATLAB session,
            % acheived by assigning return value to Constant attribute.
            je = jenv;
            if ~contains(je.Version, "Java 17")
                error("OpenRocket v23.09 requires Java 17, got '%s'", je.Version);
            end

            jcp = string(javaclasspath("-static"));
            if ~any(contains(lower(jcp), "openrocket.jar"))
                error("No 'openrocket.jar' (case-insensitive) found on static class path.")
            end

            wb = waitbar(0, "Starting OpenRocket modules..");
            gui_module = net.sf.openrocket.startup.GuiModule();
            waitbar(0.2, wb);
            plugin_module = net.sf.openrocket.plugin.PluginModule();
            waitbar(0.4, wb)
            injector = com.google.inject.Guice.createInjector([gui_module, plugin_module]);
            waitbar(0.6, wb)
            net.sf.openrocket.startup.Application.setInjector(injector);
            waitbar(0.8, wb);
            gui_module.startLoader();

            waitbar(0.9, wb);
            logger = org.slf4j.LoggerFactory.getLogger(...
                ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
            level = ch.qos.logback.classic.Level.ERROR;
            logger.setLevel(level);

            waitbar(1, wb, "Done.");
            close(wb);
            ret = true;
        end

        function dict = make_type_map
            % Create flight data type dictionary
            %   Maps name -> FlightDataType object (used by FlightDataBranch.get)
            dict = dictionary;

            import net.sf.openrocket.simulation.FlightDataType;
            % Using FlightDataType.ALL_TYPES does not always correctly recover
            % the string representations because of a localization (?) bug
            dict("Time") = FlightDataType.TYPE_TIME;
            dict("Altitude") = FlightDataType.TYPE_ALTITUDE;
            dict("Vertical velocity") = FlightDataType.TYPE_VELOCITY_Z ;
            dict("Vertical acceleration") = FlightDataType.TYPE_ACCELERATION_Z;
            dict("Total velocity") = FlightDataType.TYPE_VELOCITY_TOTAL;
            dict("Total acceleration") = FlightDataType.TYPE_ACCELERATION_TOTAL;
            dict("Position East of launch") = FlightDataType.TYPE_POSITION_X;
            dict("Position North of launch") = FlightDataType.TYPE_POSITION_Y;
            dict("Lateral distance") = FlightDataType.TYPE_POSITION_XY;
            dict("Lateral direction") = FlightDataType.TYPE_POSITION_DIRECTION;
            dict("Lateral velocity") = FlightDataType.TYPE_VELOCITY_XY;
            dict("Lateral acceleration") = FlightDataType.TYPE_ACCELERATION_XY;
            dict("Latitude") = FlightDataType.TYPE_LATITUDE;
            dict("Longitude") = FlightDataType.TYPE_LONGITUDE;
            dict("Gravitational acceleration") = FlightDataType.TYPE_GRAVITY;
            dict("Angle of attack") = FlightDataType.TYPE_AOA;
            dict("Roll rate") = FlightDataType.TYPE_ROLL_RATE;
            dict("Pitch rate") = FlightDataType.TYPE_PITCH_RATE;
            dict("Yaw rate") = FlightDataType.TYPE_YAW_RATE;
            dict("Mass") = FlightDataType.TYPE_MASS;
            dict("Motor mass") = FlightDataType.TYPE_MOTOR_MASS;
            dict("Longitudinal moment of inertia") = FlightDataType.TYPE_LONGITUDINAL_INERTIA;
            dict("Rotational moment of inertia") = FlightDataType.TYPE_ROTATIONAL_INERTIA;
            dict("CP location") = FlightDataType.TYPE_CP_LOCATION;
            dict("CG location") = FlightDataType.TYPE_CG_LOCATION;
            dict("Stability margin") = FlightDataType.TYPE_STABILITY;
            dict("Mach number") = FlightDataType.TYPE_MACH_NUMBER;
            dict("Reynolds number") = FlightDataType.TYPE_REYNOLDS_NUMBER;
            dict("Thrust") = FlightDataType.TYPE_THRUST_FORCE;
            dict("Drag force") = FlightDataType.TYPE_DRAG_FORCE;
            dict("Drag coefficient") = FlightDataType.TYPE_DRAG_COEFF;
            dict("Axial drag coefficient") = FlightDataType.TYPE_AXIAL_DRAG_COEFF;
            dict("Friction drag coefficient") = FlightDataType.TYPE_FRICTION_DRAG_COEFF;
            dict("Pressure drag coefficient") = FlightDataType.TYPE_PRESSURE_DRAG_COEFF;
            dict("Base drag coefficient") = FlightDataType.TYPE_BASE_DRAG_COEFF;
            dict("Normal force coefficient") = FlightDataType.TYPE_NORMAL_FORCE_COEFF;
            dict("Pitch moment coefficient") = FlightDataType.TYPE_PITCH_MOMENT_COEFF;
            dict("Yaw moment coefficient") = FlightDataType.TYPE_YAW_MOMENT_COEFF;
            dict("Side force coefficient") = FlightDataType.TYPE_SIDE_FORCE_COEFF;
            dict("Roll moment coefficient") = FlightDataType.TYPE_ROLL_MOMENT_COEFF;
            dict("Roll forcing coefficient") = FlightDataType.TYPE_ROLL_FORCING_COEFF;
            dict("Roll damping coefficient") = FlightDataType.TYPE_ROLL_DAMPING_COEFF;
            dict("Pitch damping coefficient") = FlightDataType.TYPE_PITCH_DAMPING_MOMENT_COEFF;
            dict("Coriolis acceleration") = FlightDataType.TYPE_CORIOLIS_ACCELERATION;
            dict("Reference length") = FlightDataType.TYPE_REFERENCE_LENGTH;
            dict("Reference area") = FlightDataType.TYPE_REFERENCE_AREA;
            dict("Vertical orientation (zenith)") = FlightDataType.TYPE_ORIENTATION_THETA;
            dict("Lateral orientation (azimuth)") = FlightDataType.TYPE_ORIENTATION_PHI;
            dict("Wind velocity") = FlightDataType.TYPE_WIND_VELOCITY;
            dict("Air temperature") = FlightDataType.TYPE_AIR_TEMPERATURE;
            dict("Air pressure") = FlightDataType.TYPE_AIR_PRESSURE;
            dict("Speed of sound") = FlightDataType.TYPE_SPEED_OF_SOUND;
            dict("Simulation time step") = FlightDataType.TYPE_TIME_STEP;
            dict("Computation time") = FlightDataType.TYPE_COMPUTATION_TIME;
        end

        function dict = make_units_map
            % Create flight units dictionary
            %   Maps name -> unit
            map = entries(openrocket.make_type_map);
            
            string_si_unit = @(fdtype) string(fdtype.getUnitGroup().getSIUnit());
            units = arrayfun(string_si_unit, map.Value);
            units(units == "°") = "deg";
            units(units == "​") = "";

            dict = dictionary(map.Key, units);
            dict("Stability margin") = "cal";
        end

        % https://www.mathworks.com/help/compiler_sdk/java/rules-for-data-conversion-between-java-and-matlab.html
        function doubles = jarr2double(jarr)
            % Convert Java ArrayList to MATLAB double array
            double_arr = javaArray("java.lang.Double", jarr.size);
            jarr.toArray(double_arr);
            doubles = double(double_arr);
            % slight speed improvment over <list = double(toArray(jarr))>;
        end

    end

    
end

