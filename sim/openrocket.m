%% OPENROCKET
% MATLAB utility class that loads an OpenRocket application and wraps
% underlying Java methods for automating OpenRocket calculations and
% simulations.
% 
% PREREQUISITES
%   - OpenRocket installation, containing jre/ and OpenRocket.jar class file
%   - Version of MATLAB <a href="https://www.mathworks.com/support/requirements/openjdk.html">compatible</a> with the Java version used by OpenRocket. 
%   At time of writing, the latest release is OR 23.09 using JDK 17, which
%   requires MATLAB 2024a or newer
%   - Set up Java class path and Java environment by running
%   openrocket_setup(...) once
%
% OVERVIEW
% The class starts the required Java modules to open an OpenRocket document,
% make modifications, and save it. The object methods are concise wrappers for
% accessing parts of the document. The class methods are concise wrappers for
% acting on parts of the document.
%
% In general,
% - All Java objects in OpenRocket are reference-type (handle), not
% copy-on-write as is typical for MATLAB. Unless an object is deep-copied,
% modifying instances will modify the underlying OpenRocket document.
% - If you are programmatically modifying the simulation, you will likely need to
% use the relevant Java object's methods directly (getOptions(), setHeight(),
% ...). They are usually self-explanatory (if somewhat cumbersome). 
% - Tab-completion works for Java object methods, but to get a complete list of
% the methods available, you can use methodsview(class(java_object)), for example:
%   or = openrocket("test.ork");
%   fins = or.component(class = "FinSet");
%   methodsview(class(fins))
%
% There are four types of OpenRocket object to pay attention to:
% - Document, which contains the Rocket, FlightConfigurations, and Simulations.
% - RocketComponent, which store the components' properties and possibly other
% components.
% - FlightConfiguration, which store settings for recovery events and motor
% mounts.
% - Simulation, which contain the settings for a given simulaation, including
% the Rocket, FlightConfiguration, and SimulationOptions.
% 
%% LISTING
% The <document> member is theoretically enough to perform any action in
% OpenRocket programmatically. These functions are defined for convenience.
%
% [Inputs and outputs in brackets are optional]
% Document methods ------------------------------------------------------------
% or = openrocket(path)
%   Load .ork file into object
% or.save()
%   Save modifications to original file
% component = or.components(Name, Value)
%   Find components by exact name (e.g. "Trapezoidal Fins", "Streamer") or
%   RocketComponent subclass (e.g. "NoseCone", "RecoveryDevice")
% config = or.get_config([identifier])
%   Get FlightConfiguration by name or number, defaulting to active one
% or.set_config(config)
%   Set active FlightConfiguration
% 
% Aerodynamic and inertial calculations ---------------------------------------
% (Unless otherwise specified, these methods use the active FlightConfiguration)
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
% 
% Class methods ---------------------------------------------------------------
% [data = ]openrocket.simulate(sim[, stop = <stop condition>][, outputs = <variable names>]);
%   Execute simulation, optionally specifying an early stop, optionally converting output data
% data = openrocket.get_data(sim[, variables])
%   retrieve data from up-to-date simulation
% openrocket.list_variables()                   
%   list possible flight data variables
% event = openrocket.get_deploy(chute, sim)
%   Get the deployment event a parachute is configured with in a specific sim
% config = openrocket.copy_config(config[, name])
%   Deep-copy a FlightConfiguration to modify sim without affecting existing configurations

% Written by Yevgeniy Gorbachev
% November 2024

%% NOTES FOR DEVELOPERS
% - The application startup is translated from <a href="https://github.com/SilentSys/orhelper/blob/master/orhelper/_orhelper.py">SilentSys/orhelper</a>, many thanks 
% for doing the annoying part.
% - There is no significant performance difference between <orhelper> and this
% class, because the bulk of the performance cost is the simulation itself.
% - Operations cannot be parallelized becuase PARFOR <a href="https://www.mathworks.com/help/parallel-computing/objects-and-handles-in-parfor-loops.html">requires</a> the data 
% passed in to support SAVE and LOAD, which Java objects in MATLAB do not
 
classdef (Sealed) openrocket < handle & matlab.mixin.Scalar
    
    properties (SetAccess = private, GetAccess = public)
        document; % OpenRocketDocument Java object
    end

    properties (Access = private)
        file;   % Java File object
        loader; % Java Loader object
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

%% STATIC UTILITIES
    methods (Static, Access = public)
        function data = simulate(sim, params)
            % Simulate an OpenRocket simulation 
            %   [data = ]openrocket.simulate(sim[, Name = Value])
            %       sim         Simulation object
            %
            %   Name-value arguments 
            %   All optional - none mutually exclusive
            %       outputs     Convert flight data - required IF [data =] output is used
            %                   String array e.g. ["Altitude", "Stability margin"]
            %                   magic value of "ALL" outputs all values
            %       stop        "apogee" to stop at apogee
            %                   positive duration to stop at specific time
            %       atmos       Custom atmospheric model
            %                   Table with columns HGT [m ASL], TMP [K], PRES [Pa]
            %       wind        Custom wind model
            %                   Table with columns HGT [m ASL], UGRD [m/s], VGRD [m/s]
            %                   U - Eastward component / V - Northward component
            %       drag        Custom drag curve
            %                   Table with columns MACH [-], DRAG [-]
            %
            %   EXAMPLES
            %       or = openrocket("data/OTIS.ork");
            %       sim = or.sims(1);
            %       % Simulate using object, up to apogee, do not output data
            %       openrocket.simulate(sim, stop = "apogee"); 
            %       Simulate and get all outputs
            %       data = openrocket.simulate(sim, output = "ALL");  

            arguments
                sim (1, 1) {openrocket.mustBeA(sim, "document.Simulation")};
                params.stop (1, :) = [];
                params.outputs (1, :) string = [];
                params.atmos table = table; 
                params.wind table = table;
                params.drag table = table;
            end

            % NOTE if this is modified, make sure the input to simulate() is
            % exactly the class Java expects (SimulationListener array) - if it
            % isn't, MATLAB will emit an inane error from a random Toolbox with
            % a simulate() function instead of an error related to the Java object.

            % import net.sf.openrocket.simulation.FlightEvent;
            listener_class = "net.sf.openrocket.simulation.listeners.SimulationListener"; 
            extension_class = "net.sf.openrocket.simulation.extension.SimulationExtension";
            
            outs = params.outputs;
            if isempty(outs) && nargout > 0
                error("Simulation output assigned, but not requested through 'outputs = ...'");
            elseif ~isempty(outs) && nargout == 0
                warning("Simulation output requested through 'outputs = ....', but not assigned.");
            end

            % Determine which extensions to add based on the optional arguments
            
            % Initialize extension array
            extensions = {};
    
            if isstring(params.stop) %&& (params.stop == "apogee")
                extensions(end+1) = OpenRocketExtensions.EarlySimulationStop(...
                    openrocket.flight_event(upper(params.stop)), 0);
            elseif isduration(params.stop)
                 extensions(end+1) = OpenRocketExtensions.EarlySimulationStop(...
                    openrocket.flight_event("IGNITION"), seconds(params.stop));
            end

            if ~isempty(params.atmos)
                openrocket.mustHaveCols(params.atmos, ["HGT", "TMP", "PRES"]);
                extensions(end+1) = OpenRocketExtensions.TabulatedAtmosphere(...
                    params.atmos.HGT, params.atmos.TMP, params.atmos.PRES);
            end

            if ~isempty(params.wind)
                openrocket.mustHaveCols(params.wind, ["HGT", "UGRD", "VGRD"]);
                extensions(end+1) = OpenRocketExtensions.TabulatedWind(...
                    params.wind.HGT, params.wind.UGRD, params.wind.VGRD);
            end

            if ~isempty(params.drag)
                openrocket.mustHaveCols(params.drag, ["MACH", "DRAG"]);
                extensions(end+1) = OpenRocketExtensions.TabulatedDrag(...
                    params.drag.MACH, params.drag.DRAG);
            end

            % convert extension cell array to many arguments to fit what List.of expects
            extensions = java.util.List.of(extensions{:});

            % Execute 
            listeners = javaArray(listener_class, 0); % incantation for Java side
            sim.copyExtensionsFrom(extensions);
            sim.simulate(listeners)

            % clear extensions to not corrupt the document
            sim.copyExtensionsFrom(java.util.List.of()); 
            
            % output data, if requested
            if ~isempty(outs) && nargout > 0
                if outs == "ALL"
                    data = openrocket.get_data(sim);
                else
                    data = openrocket.get_data(sim, outs);
                end
            end
        end

        function data = get_data(sim, variables, branch)
            % Retreive data from an OpenRocket simulation
            % data = openrocket.get_data(sim[, variables])
            %   sim         Simulation object (returned by openrocket.sims() or directly accessed)
            %   variables   (Optional) list of variables to return
            %               Defaults to add all of them
            % 
            % NOTE Does not automatically simulate outdated simulation
            % NOTE All variables are unitless or MKS
            arguments
                sim {openrocket.mustBeA(sim, "document.Simulation")};
                variables string = openrocket.types.keys;
                branch (1,1) double = 0;
            end

            mustBeMember(variables, openrocket.types.keys);
            % mustBeMember(variables, openrocket.types.keys, ...
            %     "Invalid value for <variables>. ")

            if string(sim.getStatus) ~= openrocket.status_uptodate
                error("Simulation not up to date. Run openrocket.simulate().");
            end

            variables = variables(variables ~= openrocket.time_name);
            branch = sim.getSimulatedData().getBranch(branch); 

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

        function copy = copy_config(cfg, name)
            % Deep-copy a flight configuration
            %   copy = copy_config(config[, name]);
            %       config      FlightConfiguration to copy
            %       name        (Optional) New configuration's name
            %                   Defaults to old config's name + " Copy"
            %       
            %       copy        copied FlightConfiguration

            arguments
                cfg {openrocket.mustBeA(cfg, "rocketcomponent.FlightConfiguration")}
                name (1,1) string = string(cfg.getName) + " Copy";
            end

            copy = cfg.copy([]);
            rkt = cfg.getRocket;
            configables = openrocket.search_components(rkt, class = "FlightConfigurableComponent");
            for i = 1:length(configables)
                configables(i).copyFlightConfiguration(cfg.getId, copy.getId);
            end

            if nargin == 1
                copy.setName(string(copy.getName) + " Copy");
            elseif nargin == 2
                copy.setName(name);
            end
        end

        function dev = get_deploy(chute, sim)
            % Get the deployment event for a recovery device
            %   event = or.get_deploy(chute)
            %       Use the globally selected flight configuration
            %   event = or.get_deploy(chute, sim)
            %       Use the simulation's specific flight configuration
            % Set the event's recovery configuration using the methods
            %    event.setRecoveryEvent(event.APOGEE);
            % 
            %   NOTE: Event types have to be accessed through the event object
            %   (event.APOGEE) instead of just (APOGEE) because the
            %   recoveryEvent enumeration is anonymous inside the RecoveryEvent
            %   java class
            arguments
                chute (1,1) {openrocket.mustBeA(chute, "rocketcomponent.RecoveryDevice")};
                sim (1, 1) {openrocket.mustBeA(sim, "document.Simulation")};
            end

            id = sim.getFlightConfigurationId;
            dev = chute.getDeploymentConfigurations.get(id);
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
            % or = openrocket(<path to ORK file>);
            arguments
                orkpath (1,1) string {mustBeFile};
            end

            [~, ~, ext] = fileparts(orkpath);
            if ext ~= ".ork"
                error("File has extension '%s', expected '.ork'", ext);
            end

            obj.file = java.io.File(orkpath);
            obj.loader = net.sf.openrocket.file.GeneralRocketLoader(obj.file);
            obj.document = obj.loader.load();
        end

        function save(obj)
            % Save OpenRocket document to file
            obj.saver.save(obj.file, obj.document);
        end

        function export(obj, path)
            % Save OpenRocket document to RASAero file
            %   export(ork, dest)
            %   ork     OpenRocket instance
            %   dest    Destination
            arguments
                obj (1,1) openrocket
                path (1,1) string;
            end

            [~, ~, ext] = fileparts(path);
            if ext ~= ".CDX1"
                error("RasAero export files must have extension CDX1")
            end

            file = java.io.FileOutputStream(path);
            oc = onCleanup(@() file.close());

            % These inputs are required but not used - they can stay empty
            opts = net.sf.openrocket.document.StorageOptions();
            warns = net.sf.openrocket.logging.WarningSet();
            errs = net.sf.openrocket.logging.ErrorSet();

            exporter = net.sf.openrocket.file.rasaero.export.RASAeroSaver();
            exporter.save(file, obj.document, opts, warns, errs);
        end

        function out = sims(obj, sim, nout)
            % Get simulations 
            %   sim = or.sims()       
            %       returns all simulations
            %   sim = or.sims(index)  
            %       return simulation at 1-based indices (1-based)
            %   sim = or.sims(name[, nout = length(name)])
            %       return simulations matching name(s) or <pattern> object
            %
            %       nout defaults to specify the same size, so that asking for
            %       N names expects N outputs, but can be changed for (for
            %       example) pattern-matched outputs. It will cause sims() to
            %       error on an incorrect number of outputs so that
            %       error-checking does not need to be repeatedly performed
            %       outside the funciton.
            
            arguments (Input)
                obj openrocket;
                sim (1,:) {mustBeA(sim, ["numeric", "string", "pattern"])};
                nout (1,1) double = length(sim);
            end

            sims = toArray(obj.document.getSimulations());
            if nargin == 1
                out = sims;
            elseif isnumeric(sim)
                out = sims(sim);
            elseif isstring(sim)
                names = arrayfun(@(sim) string(sim.getName()), sims);
                out = sims(matches(names, sim));
                if length(out) ~= nout
                    error("Expected %d output(s) but matched %d of%s", ...
                        nout, length(out), compose("\n\t%s", names).join(""));
                end
            else 
                error("Identifier not supported");
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

        function components = component(obj, mode, key)
            % components = obj.component(Name, Value)
            % components = obj.component(class = <class name>)
            %   class   Java class name under net.sf.openrocket.rocketcomponent
            % components = obj.component(name = <component name>)
            %   name    Component name - exact case-sensitive match
            % 
            %   EXAMPLES
            %   or = openrocket("test.ork");
            %   fins = or.component(class = "FinSet");
            %   streamer = or.component(name = "Streamer");

            arguments
                obj openrocket;
                mode (1, 1) string {mustBeMember(mode, ["name", "class"])};
                key (1, 1) string;
            end

            rkt = obj.rocket();

            switch mode
                case "name"
                    components = openrocket.search_components(rkt, "all", name = key);
                case "class"
                    components = openrocket.search_components(rkt, "all", class = key);
            end
        end

        function data = add_indicated(~, data)
            % Add flight computer indicated altitude as variable "Indicated altitude"
            % data = obj.add_indicated(data)
            %   data    (timetable)     Output with column "Air pressure"
            arguments
                ~
                data timetable {openrocket.mustHaveCols(data, "Air pressure")};
            end

            upress = data.Properties.VariableUnits{"Air pressure"};
            press = data.("Air pressure");
            if ismember("Altitude", string(data.Properties.VariableNames));
                uhgt = data.Properties.VariableUnits{"Altitude"};
            else
                uhgt = "m";
                warning("'Altitude' not presesnt or is unitless, assuming 'm'");
            end

            if upress == ""
                warning("Pressure variable does not have assigned units. Assuming Pa");
                upress = "Pa";
            end

            switch upress
                case "Pa"
                    P_mbar = press ./ 100;
                case "kPa"
                    P_mbar = press .* 10;
                case "mbar"
                    P_mbar = press;
                case "psi"
                    P_mbar = 68.94757293168361 .* press;
                otherwise
                    error("Unrecognized pressure unit '%s'", upress);
            end

            hgt_m = (1 - (P_mbar ./ 1013.25).^(0.190284)) .* 145366.45 * 0.3048;
            alt_m = hgt_m - hgt_m(1);

            % Convert output
            switch uhgt
                case "m"
                    alt = alt_m;
                case "km"
                    alt = alt_m ./ 1000;
                case "ft"
                    alt = alt_m .* 3.280839895013123;
                case "mi"
                    alt = alt_m .* 6.213711922373339e-04;
                otherwise
                    error("Unrecognized height unit '%s'", uhgt);
            end

            data.("Indicated altitude") = alt;
            data.Properties.VariableUnits{"Indicated altitude"} = uhgt;
        end

        function data = add_flutter(obj, data, modulus)
            % Add fin flutter velocity to simulation output as variable "Flutter velocity"
            % data = obj.add_flutter(data, modulus)
            %   data        (timetable)     Output with columns "Air pressure" and "Speed of sound"
            %   modulus     (double)        [Pa] fin shear modulus
            %
            % Uses flutter velocity formula from Apogee Peak of Flight 615

            arguments
                obj openrocket;
                data timetable {openrocket.mustHaveCols(data, ...
                    ["Air pressure", "Speed of sound"])}
                modulus (1,1) double {mustBePositive};
            end
            
            fins = obj.component(class = "TrapezoidFinSet");
            
            if isempty(fins)
                error("No trapezoidal fins found")
            elseif ~isscalar(fins)
                error("Multiple trapezoidal fin sets found")
            end
            % Calculations lifted from RIT-Launch-Initiative/finoptimization,
            % but applied across all times individually instead of just max. velocity
            
            % Retrieve fin dimensions
            t = fins.getThickness();
            Ls = fins.getSweep();
            Lt = fins.getTipChord();
            Lr = fins.getRootChord();
            h = fins.getHeight();

            % Constants
            P0 = 101325.3; % [Pa] Standard pressure
            k = 1.4; % [-] Specific heat ratio for air
            Cm = 24; % [-] Martins constant
            
            Cx = ((2 * Ls * Lt) + (Lt^2) + (Lr * Ls) + (Lt * Lr) + (Lr^2))/(3 * (Lt+Lr));
            eps = (Cx/Lr) - 0.25;
            Dc = (Cm * k * eps * P0)/(pi); % [Pa] "Denominator constant"
            lambda = Lt/Lr; % [-] Taper ratio
            A = h * ((Lt + Lr)/2); % [m^2] Fin area
            Ar = (h^2)/A; % [-] Aspect ratio

            P = data.("Air pressure");
            c = data.("Speed of sound");
            v_flut = c .* sqrt(modulus ./ ...
                ( ((Dc*Ar^3)/((t/Lr)^3 * (Ar + 2))) .* ((lambda + 1)/2) .* (P/P0) ));
            data.("Flutter velocity") = v_flut;


            data.Properties.VariableUnits{"Flutter velocity"} = ...
                data.Properties.VariableUnits{"Speed of sound"};
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
            arguments
                obj openrocket;
                fc {openrocket.mustBeA(fc, "aerodynamics.FlightConditions")};
            end

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
            arguments
                obj openrocket;
                fc {openrocket.mustBeA(fc, "aerodynamics.FlightConditions")};
            end

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
            arguments
                obj openrocket;
                state (1,1) string {mustBeMember(state, ["LAUNCH", "BURNOUT"])};
            end

            cfg = obj.get_config();
            switch state
                case "LAUNCH"
                    data = openrocket.masscalc.calculateLaunch(cfg);
                case "BURNOUT"
                    data = openrocket.masscalc.calculateBurnout(cfg);
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
            arguments
                obj openrocket;
                state (1,1) string {mustBeMember(state, ["LAUNCH", "BURNOUT"])};
                fc {openrocket.mustBeA(fc, "aerodynamics.FlightConditions")};
            end

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

            % wb = waitbar(0, "Starting OpenRocket modules..");
            gui_module = net.sf.openrocket.startup.GuiModule();
            % waitbar(0.2, wb);
            plugin_module = net.sf.openrocket.plugin.PluginModule();
            % waitbar(0.4, wb)
            injector = com.google.inject.Guice.createInjector([gui_module, plugin_module]);
            % waitbar(0.6, wb)
            net.sf.openrocket.startup.Application.setInjector(injector);
            % waitbar(0.8, wb);
            gui_module.startLoader();

            % waitbar(0.9, wb);
            logger = org.slf4j.LoggerFactory.getLogger(...
                ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
            level = ch.qos.logback.classic.Level.ERROR;
            logger.setLevel(level);

            % waitbar(1, wb, "Done.");
            % close(wb);
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

        function ev = flight_event(name)
            arguments
                name (1,1) string;
            end

            % internal methods are converted to $<name> instead of .<name> so
            % we have to do this garbage
            ev = javaMethod("valueOf", "net.sf.openrocket.simulation.FlightEvent$Type", name);
        end

        % https://www.mathworks.com/help/compiler_sdk/java/rules-for-data-conversion-between-java-and-matlab.html
        
        function doubles = jarr2double(jarr)
            % Convert Java ArrayList to MATLAB double array
            double_arr = javaArray("java.lang.Double", jarr.size);
            jarr.toArray(double_arr);
            doubles = double(double_arr);
            % slight speed improvment over <list = double(toArray(jarr))>;
        end

        function tf = mustBeA(input, classname, prefix)
            % Validate that the input is or inherits from the specified class name 
            % Behaves like <mustBeA>, but prepends net.sf.openrocket 
            % not necessary for MATLAB types, which support direct validation
            arguments
                input
                classname (1,1) string;
                prefix (1,1) string = "";
            end
            fullname = "net.sf.openrocket." + classname;
            tf = true;
            if ~isa(input, fullname)
                ME = MException("openrocket:invalid_type", ...
                    sprintf("\n%sExpected %s\nGot %s", prefix, fullname, class(input)));
                throwAsCaller(ME);
            end
        end

        function mustHaveCols(input, cols)
            arguments
                input;
                cols (1,:) string;
            end

            err_id = "openrocket:invalidTableInput";
            notpresent = setdiff(cols, input.Properties.VariableNames);
            if ~isempty(notpresent)
                mex = MException(err_id, "Required table columns %s not present", ...
                    mat2str(notpresent));
                throwAsCaller(mex);
            end
        end

        function components = search_components(parent, depth, mode, key)
            % Search rocket for components
            % components = search_components(parent, depth, Name = Value)
            %   parent  Rocket object (or subcomponent)
            %   depth   "current" or "all" for immediate children or all children
            % components = search_components(..., class = <class name>)
            %   class   Java class name under net.sf.openrocket.rocketcomponent
            % components = search_components(..., class = <class name>)
            %   name    Component name - exact case-sensitive match
            % 
            %   EXAMPLES
            %   or = openrocket("test.ork");
            %   fins = openrocket.search_components(or.rocket(), "all", class = "FinSet");
            %   streamer = openrocket.search_components(or.rocket(), "all", name = "Streamer");

            arguments
                parent {openrocket.mustBeA(parent, ...
                    "rocketcomponent.RocketComponent")};
                depth (1, 1) string {mustBeMember(depth, ["current", "all"])};
                mode (1, 1) string {mustBeMember(mode, ["name", "class"])};
                key (1, 1) string;
            end 

            
            switch depth
                case "current"
                    components = parent.getChildren().toArray();
                case "all"
                    components = parent.getAllChildren().toArray();
            end

            switch mode
                case "name"
                    match = string(components) == key;
                case "class"
                    java_class_name = "net.sf.openrocket.rocketcomponent." + key;
                    match = false(size(components));
                    % for-loop required because isa() does not vectorize over an Array(List)
                    for i_part = 1:numel(components)
                        match(i_part) = isa(components(i_part), java_class_name);
                    end
            end

            components = components(match);
        end
    end

end

