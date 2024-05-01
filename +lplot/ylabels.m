% Add primary and secondary Y-axis labels to a plot
% ylabels(primary, unit)
%   primary:    axis label, passed to ylabel
%   unit:       unit label -- MATLAB figures out exponent and appends unit
function ylabels(primary, unit)
    arguments
        primary (1,1) string;
        unit (1,1) string;
    end
    import lplot.utl.get_ax_text;

    yaxis = gca().YAxis;
    unit = get_ax_text(yaxis, unit);
    
    ylabel(primary);
    ysecondarylabel(unit);
end
