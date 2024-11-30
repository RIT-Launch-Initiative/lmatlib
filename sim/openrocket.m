%% OPENROCKET
% MATLAB implementation of the <a
% href="https://github.com/SilentSys/orhelper">orhelper</a> Python script
% that provides access to OpenRocket documents for easier simulation
% automation and control system design.
% 
% PREREQUISITES
%   - OpenRocket installation, containing jre/ and OpenRocket.jar class file
%   - Version of MATLAB <a href="https://www.mathworks.com/support/requirements/openjdk.html">compatible</a> with the Java version used by OpenRocket
%   (Current OR 23.09 uses JDK 17, which requires MATLAB 2024a or newer)
%   - Set up Java class path and Java environment using openrocket_setup(<path to openrocket>)
%   
% NOTES 
%   - The OpenRocket class file has several components that duplicate classes
%   MATLAB uses, so an SL4J error will appear on startup. This appears to
%   practically not matter.
%   - Typically, OpenRocket will be located at C:\Program Files\OpenRocket
%   - Many of OpenRocket's features are, in one way or another, accessible
%   through the <document> property of this class. 

% Written by Yevgeniy Gorbachev
% November 2024
classdef openrocket < handle
    
    properties (SetAccess = protected, GetAccess = public)
        document; % OpenRocketDocument Java object
    end

    properties (Constant, Access = public)
        types = openrocket.make_type_map();
        units = openrocket.make_units_map();
    end

    properties (Constant, Access = protected)
        rename = openrocket.make_rename_map();
        started = openrocket.start(); % dummy constant to ensure start() is called once
        saver = net.sf.openrocket.file.GeneralRocketSaver();

        % Magic values
        status_uptodate = "UPTODATE"; % indicates simulation is up to date
        listener_class = "net.sf.openrocket.simulation.listeners.SimulationListener"; 
        % ^ data type for SimulationListener to make an empty array
        time_name = "Time"; % Name of the time column
        reco_name = openrocket.translate_event("Recovery device deployment"); 
        % ^ Name of a recovery device deployment
    end

    properties (Access = protected)
        file;   % Java File object
        loader; % Java Loader object
    end


    methods (Static, Access = public)
        function simulate(sim)
            % Simulate an OpenRocket simulation with no listeners attached
            %   sim     Simulation object
            if string(sim.getStatus) ~= openrocket.status_uptodate % So that repeated calls don't waste time
                sim.simulate(javaArray(openrocket.listener_class, 0)); % NOTE simulates with no listeners
            end
        end

        % https://www.mathworks.com/help/compiler_sdk/java/rules-for-data-conversion-between-java-and-matlab.html
        function data = get_data(sim, variables)
            % Retreive data from an OpenRocket simulation
            % data = get_data(sim, variables)
            %   sim         Simulation object (returned by openrocket.sims() or directly accessed)
            %   variables   (Optional) list of variables to return
            %               Defaults to add all of them
            % 
            % NOTE All variables are unitless or MKS
            arguments
                sim
                variables string = openrocket.types.keys;
            end

            variables = variables(variables ~= openrocket.time_name);
            openrocket.simulate(sim); % Ensure simulation is up to date
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
            names = arrayfun(@(ev) openrocket.translate_event(ev.getType), evs);
            times = arrayfun(@(ev) ev.getTime, evs);
            reco_idx = find(names == openrocket.reco_name);
            if length(reco_idx) == 2
                names(reco_idx(1)) = "DROGUE";
                names(reco_idx(2)) = "MAIN";
            end
            data.Properties.Events = eventtable(seconds(times), EventLabels = names); 
            data.Properties.VariableUnits = openrocket.units(data.Properties.VariableNames);
            data.Properties.VariableNames = openrocket.rename(data.Properties.VariableNames);
            data.Properties.VariableContinuity = repmat("continuous", 1, width(data));
            data.Properties.Description = "openrocket"; % Identify as OpenRocket import
        end

        function vars = list_variables
            % List available flight data types
            vars = openrocket.types.keys;
        end
    end

    methods (Static, Access = protected)
        function ret = start()
            % Create barebones OpenRocket application. This is intended to be
            % called exactly once in a MATLAB session, acheived by assigning return
            % value to Constant attribute.
            je = jenv;
            if ~contains(je.Version, "Java 17")
                error("OpenRocket v23.09 requires Java 17, got '%s'", je.Version);
            end

            jcp = string(javaclasspath("-static"));
            if ~any(contains(lower(jcp), "openrocket.jar"))
                error("No 'openrocket.jar' (case-insensitive) found on static class path.")
            end

            gui_module = net.sf.openrocket.startup.GuiModule();
            plugin_module = net.sf.openrocket.plugin.PluginModule();

            injector = com.google.inject.Guice.createInjector([gui_module, plugin_module]);
            net.sf.openrocket.startup.Application.setInjector(injector);
            gui_module.startLoader();

            logger = org.slf4j.LoggerFactory.getLogger(...
                ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
            level = ch.qos.logback.classic.Level.ERROR;
            logger.setLevel(level);
            ret = true;
        end

        function dict = make_type_map
            % Create flight data type dictionary
            types = net.sf.openrocket.simulation.FlightDataType.ALL_TYPES;
            keys = string(types);
            dict = dictionary(keys, types);
        end

        function dict = make_units_map
            types = net.sf.openrocket.simulation.FlightDataType.ALL_TYPES;
            keys = string(types);
            units = arrayfun(@(fdt) string(fdt.getUnitGroup().getSIUnit()), types);

            units(units == "°") = "deg";
            units(units == "​") = "";

            dict = dictionary(keys, units);
            dict("Stability margin calibers") = "cal";
        end

        function dict = make_rename_map
            dict = dictionary(openrocket.types.keys, openrocket.types.keys);
            dict("Stability margin calibers") = "Stability margin";
        end

        function doubles = jarr2double(jarr)
            % Convert Java ArrayList to MATLAB double array
            double_arr = javaArray("java.lang.Double", jarr.size);
            jarr.toArray(double_arr);
            doubles = double(double_arr);
            % slight speed improvment over <list = double(toArray(jarr))>;
        end

        function str = translate_event(ev)
            str = upper(string(ev));
            str = strrep(str, " ", "_");
            str = strrep(str, "-", "_");
        end
    end
    
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
            %   sims()       returns all simulations
            %   sims(index)  return simulation at index (1-based)
            %   sims(name)   return all simulations matching name

            sims = toArray(obj.document.getSimulations());
            if nargin == 1
                out = sims;
            elseif isnumeric(sim)
                out = sims(sim);
            elseif isstring(sim)
                names = arrayfun(@(sim) string(sim.getName()), sims);
                out = sims(names == sim);
            end
        end

        function rkt = rocket(obj)
            % Get Rocket object from document
            rkt = obj.document.getRocket();
        end

        function comp = component(obj, name)
            % Search entire Rocket for component name
            %   component = or.component(name)

            parts = obj.rocket().getAllChildren();
            parts = toArray(parts);
            names = string(parts);
            found = find(names == name);
            if isempty(found)
                error("Could not find component '%s' in rocket", name);
            elseif length(found) >= 2
                error("Found %d components matching '%s'", length(found), name);
            end
            comp = parts(found);
        end
    end
end

