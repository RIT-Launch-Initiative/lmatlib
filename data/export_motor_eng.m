%% Write .eng motor definition file
% export_motor_eng(path, data)
% Inputs
%   path    (string)    destination path
%   data    (timetable) thrust history and metadsata
%                       .Properties.Description is written as a comment in the beginning
%                       .Properties.UserData is a structure like the outupt of import_motor_eng

function export_motor_eng(path, data)
    arguments (Input)
        path (1,1) string;
        data (:,1) timetable;
    end

    required_col = "Thrust";
    required_unit = "N";
    required_fields = ["name", "diameter_mm", "length_mm", "delays_sec", ...
        "propmass_kg", "wetmass_kg", "manufacturer"];

    names = data.Properties.VariableNames;
    units = data.Properties.VariableUnits;
    
    defs = data.Properties.UserData;
    fields = string(fieldnames(defs));
    notpresent = setdiff(required_fields, fields);

    % Validate input
    assert(names{1} == required_col, ...
        "Thrust data must have one column named '%s'", required_col);
    assert(isscalar(units) && units{1} == required_unit, ...
        "Invalid input: Thrust data must have VariableUnits of '%s'", required_unit);
    assert(all(data.Thrust >= 0), "All thrust points must be non-negative");
    assert(seconds(data.Time(1)) > 0, "Time must start after 0 seconds");
    assert(data.Thrust(end) == 0, "Last thrust point must be exactly 0");
    assert(isempty(notpresent), "Required fields of UserData not present: %s", ...
        mat2str(notpresent));

    header_items = repmat("", 1, length(required_fields));
    for i_item = 1:length(required_fields)
        header_items(i_item) = defs.(required_fields(i_item));
    end
    header_line = header_items.join(" ");

    fid = fopen(path, "w");
    fprintf(fid, "%s\n", header_line);
    fprintf(fid, compose(" %f %f", seconds(data.Time), data.Thrust).join(char(10)));
    fclose(fid);
end
