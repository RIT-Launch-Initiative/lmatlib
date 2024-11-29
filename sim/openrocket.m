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
%   - Set up Java class path and Java environment using openrocket.setup(...)
%   
% NOTES 
%   - Typically, OpenRocket will be located at C:\Program Files\OpenRocket, and
%   that is the default used by setup()
classdef openrocket < handle
    
    % NOTES for DEVELOPERS
    % - Avoid using "import" to abbreviate Java package access 
    %   - it is not consistently applicable because some invocations like
    %   javaMethod(..) and javaArray(..) require the full name anyway
    %   - All import statements are loaded before class methods can be used, so
    %   a first-time user would not be able to run the setup() function 

    properties (SetAccess = protected, GetAccess = public)
        started (1,1) logical = false; % OR started?
        ork (1,1) string = ""; % Path to .ork file
        loaded (1,1) logical = false; % .ork loaded?
    end

    properties (Access = public)
    % properties (Access = protected)
    % Eventually will protect these, but exposed to user for now
        document;
        loader;
        saver;
    end

    methods (Static, Access = public)
        function setup(openrocket_path)
            % setup(openrocket_path) - folder of OpenRocket installation 
            %   Contains OpenRocket.jar and jre/
            %   Defaults to C:\Program Files\OpenRocket
            % NOTE: Modifies $prefdir/javaclasspath.txt
            % NOTE: Modifies Java environment - ensure MATLAB supports OR's Java runtime
            arguments
                openrocket_path (1,1) string = "C:\Program Files\OpenRocket";
            end

            if ~isfolder(openrocket_path)
                error("Directory '%s' not found", java_path);
            end
            jenv_path = fullfile(openrocket_path, "jre");
            if ~isfolder(jenv_path)
                error("'jre' not found under %s", openrocket_path);
            end
            jar_path = fullfile(openrocket_path, "OpenRocket.jar");
            if ~isfile(jar_path)
                error("'OpenRocket.jar' not found under %s\n", openrocket_path);
            end
            
            jenv(jenv_path);
            jcp = fullfile(prefdir, "javaclasspath.txt");
            fid = fopen(jcp, "w");
            % Write <before> to pre-pend contents of javaclasspath to MATLAB static path
            % This prevents a conflicting (old) version of Guice from breaking the startup sequence
            % https://stackoverflow.com/questions/16366059/best-way-to-override-matlabs-default-static-javaclasspath
            fwrite(fid, "<before>");
            fwrite(fid, jar_path);
            fclose(fid);

            jenv
            fprintf("OpenRocket class path written to %s\n", jcp);
            fprintf("Restart MATLAB to apply changes\n");
        end

    end
    
    methods (Access = public)
        function obj = openrocket(orkpath)
            % Construct openrocket object with specified file
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

            obj.ork = orkpath;
        end

        function start(obj)
            % Start OpenRocket
            gui_module = net.sf.openrocket.startup.GuiModule();
            plugin_module = net.sf.openrocket.plugin.PluginModule();

            injector = com.google.inject.Guice.createInjector([gui_module, plugin_module]);
            net.sf.openrocket.startup.Application.setInjector(injector);
            gui_module.startLoader();

            logger = org.slf4j.LoggerFactory.getLogger(...
                ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
            level = ch.qos.logback.classic.Level.ERROR;
            logger.setLevel(level);

            java_file = java.io.File(obj.ork);
            obj.loader = net.sf.openrocket.file.GeneralRocketLoader(java_file);
            obj.saver = net.sf.openrocket.file.GeneralRocketSaver();

            obj.started = true;
        end

        function load(obj)
            % Load OpenRocket document
            if ~obj.started
                error("OpenRocket not started");
            end
            obj.document = obj.loader.load();
        end

        function save(obj)
            % Save OpenRocket document
            if ~obj.started
                error("OpenRocket not started");
            end
            java_file = java.io.File(obj.ork);
            obj.saver.save(java_file, obj.document);
        end
    end
end
