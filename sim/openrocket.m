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
classdef openrocket < handle
    
    properties (SetAccess = protected, GetAccess = public)
        document; % OpenRocketDocument Java object
    end

    properties (Constant, Access = protected)
        saver = net.sf.openrocket.file.GeneralRocketSaver();
        types = dictionary(...
            string(net.sf.openrocket.simulation.FlightDataType.ALL_TYPES), ...
            net.sf.openrocket.simulation.FlightDataType.ALL_TYPES);
    end

    properties (Access = protected)
        path (1,1) string = ""; % Path to OpenRocket file
        loader; % Loader object
    end


    methods (Static, Access = public)
        function simulate(sim)
            % Simulate an OpenRocket simulation with no listeners attached
            %   sim     Simulation object
            if string(sim.getStatus) ~= "UPTODATE"
                listener_dt = "net.sf.openrocket.simulation.listeners.SimulationListener";
                sim.simulate(javaArray(listener_dt, 0)); % NOTE simulates with no listeners
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
            % NOTE Does not assign units, export events, or identify table as OR import for plot_openrocket()
            arguments
                sim
                variables string = openrocket.types.keys;
            end

            variables = variables(variables ~= "Time");
            openrocket.simulate(sim);
            branch = sim.getSimulatedData().getBranch(0); % NOTE assumes branch 0

            time = openrocket.jarr2double(branch.get(openrocket.types("Time")));
            data = timetable(Size = [length(time) length(variables)], ...
                VariableTypes = repmat("double", 1, length(variables)), ...
                VariableNames = variables, ...
                RowTimes = seconds(time)); % NOTE assumes time is in seconds

            for idx = 1:length(variables)
                col = variables(idx);
                jarr = branch.get(openrocket.types(col));
                if isempty(jarr)
                    data = removevars(data, col);
                else
                    data.(col) = openrocket.jarr2double(branch.get(openrocket.types(col)));
                end
            end
        end

        function vars = list_variables
            vars = openrocket.types.keys;
        end
    end

    methods (Static, Access = protected)
        function doubles = jarr2double(jarr)
            double_arr = javaArray("java.lang.Double", jarr.size);
            jarr.toArray(double_arr);
            doubles = double(double_arr);
            % slight speed improvment over <list = double(toArray(jarr))>;
        end
    end
    
    methods 
        function obj = openrocket(orkpath)
            % Construct openrocket object and load program
            arguments
                orkpath (1,1) string
            end

            [~, ~, ext] = fileparts(orkpath);
            if ~isfile(orkpath)
                error("File '%s' does not exist", orkpath);
            end
            if ext ~= ".ork"
                error("File has extension '%s' instead of '.ork'", ext);
            end

            je = jenv;
            if ~contains(je.Version, "Java 17")
                error("OpenRocket v23.09 requires Java 17, got '%s'", je.Version);
            end

            jcp = string(javaclasspath("-static"));
            if ~any(contains(lower(jcp), "openrocket.jar"))
                error("No 'openrocket.jar' (case-insensitive) found on static class path.")
            end

            obj.path = orkpath;

            gui_module = net.sf.openrocket.startup.GuiModule();
            plugin_module = net.sf.openrocket.plugin.PluginModule();

            injector = com.google.inject.Guice.createInjector([gui_module, plugin_module]);
            net.sf.openrocket.startup.Application.setInjector(injector);
            gui_module.startLoader();

            logger = org.slf4j.LoggerFactory.getLogger(...
                ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
            level = ch.qos.logback.classic.Level.ERROR;
            logger.setLevel(level);

            java_file = java.io.File(obj.path);
            obj.loader = net.sf.openrocket.file.GeneralRocketLoader(java_file);
            obj.document = obj.loader.load();
        end

        function save(obj)
            % Save OpenRocket document
            java_file = java.io.File(obj.path);
            obj.saver.save(java_file, obj.document);
        end

        function out = sims(obj, sim)
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
    end
end

