function ph = timeplot(timetab, name, varargin)
    arguments
        timetab timetable;
        name (1,1) string;
    end
    arguments (Repeating)
        varargin;
    end

    import lplot.ylabels;

    names = timetab.Properties.VariableNames;
    units = timetab.Properties.VariableUnits;
    tt_index = find(names == name, 1, "first");
    if isempty(tt_index)
        error("Name %s not present in timetable", name);
    end

    ph = plot(timetab.Time, timetab.(name), varargin{:});

    if isempty(units) || isempty(units(tt_index))
        ylabel(name);
    else
        ylabels(name, units(tt_index));
    end
end
