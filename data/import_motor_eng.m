%% Read .eng motor definition file
% [data] = import_motor_eng(path)
% Inputs
%   path    (string)    Path to .eng file
% Outputs
%   data    (timetable) Thrust curve with columns ("Time", "Thrust"), includes associated metadata:
%                       .Properties.Description includes the comments and header line
%                       .Properties.UserData is a structure:
%                           .name       (string) motor name
%                           .diameter_mm (double) motor diameter
%                           .length_mm  (double) motor length
%                           .delays_sec (string) delay numbers separated by '-', or 'P' for plugged
%                           .propmass_kg (double) propellant mass
%                           .wetmass_kg (double) total mass
%                           .manufacturer (string) NAR manufacturer abbreviation

function [data] = import_motor_eng(path)
    arguments (Input)
        path (1,1) string {mustBeFile};
    end
    arguments (Output)
        data (:,1) timetable;
    end
    
    % readtable() is fine on its own for an ad-hoc script but anything widely
    % reusable should explicitly specify the format using
    % delimitedTextImportOptions
    expected_ext = ".eng";
    expected_num_defs = 7;


    [~, ~, ext] = fileparts(path);
    if ext ~= expected_ext
        warning("Expected extension '%s', got '%s'", expected_ext, ext);
    end
    
    lines = readlines(path);
    isblank = startsWith(lines, ";") | (lines == "");
    i_header = find(~isblank, 1, "first");
    if i_header == length(lines)
        error("Malformed motor definition: No data after header line");
    end

    defn_line = lines(i_header);


    defn_fields = defn_line.split(" ");
    if length(defn_fields) ~= expected_num_defs
        error("Malformed motor definition file: " + ...
            "Definition line does not have the %d required fields ""%s""", ...
            expected_num_defs, defn_line);
    end

    motordefs.name = defn_fields(1);
    motordefs.diameter_mm = str2double(defn_fields(2));
    motordefs.length_mm = str2double(defn_fields(3));
    motordefs.delays_sec = defn_fields(4);
    motordefs.propmass_kg = str2double(defn_fields(5));
    motordefs.wetmass_kg = str2double(defn_fields(6));
    motordefs.manufacturer = defn_fields(7);


    opts = delimitedTextImportOptions;
    opts.CommentStyle = ";";
    opts.Delimiter = " ";
    opts.LeadingDelimitersRule = "ignore";
    opts.TrailingDelimitersRule = "ignore";
    opts.EmptyLineRule = "skip";
    opts.ExtraColumnsRule = "error";
    opts.VariableNames = ["Time", "Thrust"];
    opts.DataLines = i_header+1;
    opts.VariableTypes = ["double", "double"];

    data = readtable(path, opts);
    data.Time = seconds(data.Time);
    data = table2timetable(data);

    if ~issorted(data.Time, 1, "ascend");
        error("Malformed motor definition: Time is not ascending")
    end

    if data.Thrust(end) ~= 0
        error("Malformed motor definition: last thrust point must be zero");
    end

    data.Properties.VariableUnits = "N"; % standardized
    data.Properties.UserData = motordefs;
    % add comments if there are any
    if i_header > 1
        commentlines = strip(lines(1:(i_header-1)), ";");
        data.Properties.Description = join(commentlines, char(10)); % join by newline
    end
end
