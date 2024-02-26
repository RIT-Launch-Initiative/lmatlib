% Plot one or two variables between given event labels
% plot_interval(table, first, second,  ... 
%   start_ev = "EVENT", end_ev = "EVENT" [, markers = ["EVENT1", "EVENT2", ...]])
%   table:      Table imported from OpenRocket
%   first:      variable to plot on the left axis
%   second:     (optional) variable to plot on the right axis
%   start_ev:   starting event label (must occur exactly once in the table)
%   end_ev:     ending event label (must occur exactly once in the table)
%   markers:    (optional) additonal markers between the start and end to include
%               The start and end markers are always present
%               Markers beyond the bounds start/end_ev are ignored
% 
% Example: plot_interval(omen, "Altitude", ...
%               start_ev = "LAUNCH", end_ev = "GROUND_HIT", labels = ["BURNOUT", "DROGUE", "MAIN"]);
% Example: plot_interval(omen, "Angle of attack", "Stability margin", ...
%               start_ev = "LAUNCHROD", end_ev = "APOGEE", labels = "BURNOUT");
function [pl, pr] = plot_interval(table, first, second, interval, markers)
    arguments
        table timetable;
        first (1,1) string;
        second (1,1) string = "";
        interval.start_ev (1,1) string;
        interval.end_ev (1,1) string;
        markers.labels (1, :) string = [];
    end

    plot_range = timerange(eventfilter(interval.start_ev), ...
        eventfilter(interval.end_ev), "closed");
    table = table(plot_range, :)

    % get selected columns
    if second == ""
        table = table(:, first);
    else
        table = table(:, [first, second]);
        yyaxis left;
    end
    
    % plot left axis
    plot(table, first);
    ylabel(first);
    add_unit("y", table.Properties.VariableUnits{1});

    % plot right axis if present
    if second ~= ""
        yyaxis right;

        plot(table, second);
        ylabel(second);
        add_unit("y", table.Properties.VariableUnits{2});

        yyaxis left;
    end

    xlabel("Time");

    % get markers that are selected and in range
    plot_time_range = timerange(table.Time(1), table.Time(end), "closed");
    ev_in_range = table.Properties.Events(plot_time_range, :);
    selected_ev = ismember(ev_in_range.EventLabels, ...
        [markers.labels, interval.start_ev, interval.end_ev]);

    % plot markers
    xl = xline(ev_in_range.Time(selected_ev), "-k", ...
        ev_in_range.EventLabels(selected_ev), ...
        Interpreter = "none");
    xl(end).LabelHorizontalAlignment = "left";
end

