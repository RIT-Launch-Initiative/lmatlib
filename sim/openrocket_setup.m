% openrocket_setup(path) - set up OpenRocket-MATLAB connection (run only once ever)
%   <path> is the absolute path to the OpenRocket installation 
%   Contains OpenRocket.jar and jre/
%   
% Example: openrocket_setup("C:/Program Files/OpenRocket")
% 
% NOTE: Overwrites $prefdir/javaclasspath.txt
% NOTE: Modifies Java environment - ensure MATLAB <a href="https://www.mathworks.com/support/requirements/openjdk.html">supports</a> OR's Java runtime
% NOTE: Some of OpenRocket's classes override MATLAB classes, which may cause
% errors (but no breaking changes are known of at time of writing)
function openrocket_setup(openrocket_path)
    arguments
        openrocket_path (1,1) string;
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

