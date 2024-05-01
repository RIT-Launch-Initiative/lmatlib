% Import the CSV file output by OpenRocket into a MATLAB timetable
% openrocket_table = import_or(path, verbose = false)
%   path:       path to CSV, txt, or other file readable by readtable()
%   verbose:    print recognzied variables/events/etc.
%               false by default
% 
% Export notes:
% - Include headers 
% - Include flight events
% - Use "," as the field separator
% - Use "#" as the comment character

function timetab = import_or(path, opts)
    arguments
        path (1,1) string;
        opts.verbose (1,1) logical = false;
        % structure in the arguments block gives its fields as the names for name=value arguments

        % TODO extra opts to specify comment and field-separation character instead of hard-coding
    end
    lines = readlines(path);

    assert(startsWith(lines(1), "#"), "Incorrect comment character: expected #");

    %% Create header
    header_sw = "# Time";
    header_fmt = "([A-Z].*) \((.*)\)"; % column header format (regex)
    % captures "<Name> (<unit>)"

    head_line_n = find(startsWith(lines, header_sw), 1, "first");
    if opts.verbose 
        fprintf("Header found at line %d\n", head_line_n);
    end
    head_strings = split(lines(head_line_n), ",");
    head_matches = regexp(head_strings, header_fmt, "tokens");

    % unpack cell array into (row) name and unit string arrays
    names = strings(length(head_matches), 1);
    units = strings(length(head_matches), 1);
    for i = 1:length(head_matches)
        if isempty(head_matches{i}) % doesn't follow field name format, corrupted file
                error("Could not match column header %d:\n\t""%s""", i, head_strings(i));
        end
        
        pair = head_matches{i}{1};
        if pair(1) == "Stability margin calibers"
            pair(1) = "Stability margin";
            pair(2) = "cal.";
        end

        names(i) = pair(1);
        units(i) = pair(2);
    end

    % status output
    if opts.verbose 
        disp_table = table(names, units);
        fprintf("Recognized %d columns and units\n\n", length(names));
        disp(disp_table);
    end
    
    %% Get data
    tab = readtable(path, CommentStyle = "#");
    if opts.verbose 
        fprintf("Read %d rows of data.\n", height(tab));
    end

    tab.Properties.VariableNames = names;
    tab.Properties.VariableUnits = units;

    % convert numbers to times
    for i = find(units == "s")'
        tab.(names(i)) = seconds(tab.(names(i)));
    end
    for i = find(units == "min")'
        tab.(names(i)) = seconds(60 * tab.(names(i)));
    end

    timetab = table2timetable(tab);

    %% Get event times
    event_sw = "# Event";
    event_fmt = "# Event ([A-Z_]+) occurred at t=([\d.]+) seconds"; % event line format (regex)
    % captures "Event <EVENT_NAME> occurred at t=<1234.56> seconds"

    event_strings = lines(startsWith(lines, event_sw));
    event_matches = regexp(event_strings, event_fmt, "tokens");
    event_count = length(event_matches);

    if isempty(event_strings)
        error("No events detected");
    end
    
    % unpack cell aray from event_matches into name and time arrays
    ev_names = strings(event_count, 1); 
    ev_times = duration(NaN(event_count, 3));
    for i = 1:event_count
        if isempty(event_matches{i})
            warning("Failed to match name and time in event %d: ""%s""", ...
                i, event_strings(i));
        else 
            pair = event_matches{i}{1};
            ev_names(i) = pair(1);

            % find closest time to t=...
            % because of roundoff (or other things), event time might not exactly match table time
            ev_time_raw = seconds(str2double(pair(2)));
            [~, ev_idx] = min(abs(timetab.Time - ev_time_raw)); % closest time
            ev_times(i) = timetab.Time(ev_idx);
        end
    end


    %% Name recovery devices
    reco_evname = "RECOVERY_DEVICE_DEPLOYMENT";

    is_reco = ev_names == reco_evname;
    switch sum(is_reco)
        case 0
            warning("Did not find %s", reco_evname);
        case 1
            ev_names(is_reco) = "MAIN";
            if opts.verbose
                fprintf("Converted one %s to MAIN at %.2f seconds\n", ...
                    reco_evname, seconds(ev_times(is_reco)));
            end
        case 2
            reco_name_order = ["DROGUE"; "MAIN"];
            ev_names(is_reco) = reco_name_order;
            reco_times = seconds(ev_times(is_reco));
            if opts.verbose
                fprintf("Converted first %s to %s at %.2f seconds\n", ...
                    reco_evname, reco_name_order(1), reco_times(1));
                fprintf("Converted second %s to %s at %.2f seconds\n", ...
                    reco_evname, reco_name_order(2), reco_times(2));
            end
        otherwise
            ats = compose("\tat %.2f seconds", seconds(ev_times(is_reco)));
            warning("Cannot disambiguate %d deployment events\n%s", ...
                sum(is_reco), join(ats, newline));
    end

    if numel(ev_names) ~= numel(unique(ev_names))
        warning("Non-unique event names cannot be used to delimit plot bounds");
    end

    event_table = eventtable(ev_times, EventLabels = ev_names);
    if opts.verbose 
        fprintf("Recognized %d events\n\n", height(event_table));
        disp(event_table);
    end


    % Identify this as an OpenRocket import
    timetab.Properties.Events = event_table;
    timetab = addprop(timetab, "is_openrocket_import", "table");
    timetab.Properties.CustomProperties.is_openrocket_import = true;
end
