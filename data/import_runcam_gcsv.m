function data = import_runcam_gcsv(path)
    %% data = import_runcam_gcsv(path)
    %   path    string      path to .gcsv file
    %   
    %   data    timetable   Timetable with columns "r(x/y/z)", "a(x/y/z)"
    %                       Properties.VariableUnits as "rad/s" or "m/s^2"
    % Supports files containing gyroscope and accelerometer (not magnetometer) data
    % See the <a href="https://docs.gyroflow.xyz/app/technical-details/gcsv-format">Gyroflow documentation</a>
    arguments
        path (1,1) string {mustBeFile};
    end

    %% MAGIC VALUES
    nheader = 15; % How many metadata lines 
    nvars = 1 + 3 + 3; % How many data columns
    numeric_fields = ["tscale", "gscale", "ascale"]; % Numeric metadata

    dataopts = detectImportOptions(path, ...
        FileType = "delimitedtext", ...
        NumHeaderLines = nheader, ...
        ReadVariableNames = false, ...
        ExpectedNumVariables = nvars, ...
        VariableNamesLine = nheader+1);

    headeropts = delimitedTextImportOptions(NumVariables = 2, ...
        ExtraColumnsRule = "ignore", ...
        VariableTypes = ["string", "string"], ...
        VariableNames = ["field", "value"], ...
        Delimiter = dataopts.Delimiter, ...
        DataLines = [2 nheader]);

    info = readtable(path, headeropts);
    interleaved = [info.field'; info.value'];
    metadata = struct(interleaved{:});

    for field = numeric_fields
        metadata.(field) = str2double(metadata.(field));
    end

    rawdata = readtable(path, dataopts);
    times = seconds(rawdata.t * metadata.tscale);
    gdata = rawdata{:, ["rx", "ry", "rz"]} * metadata.gscale;
    adata = rawdata{:, ["ax", "ay", "az"]} * metadata.ascale * 9.807;

    data = array2timetable([gdata adata], RowTimes = times, ...
        VariableNames = ["rx", "ry", "rz", "ax", "ay", "az"]);
    data.Properties.VariableUnits = ["rad/s", "rad/s", "rad/s", "m/s^2", "m/s^2", "m/s^2"];
    data.Properties.UserData = metadata;
end
