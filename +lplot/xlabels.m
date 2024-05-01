% Add primary and secondary X-axis labels to a plot
% xlabels(primary, unit)
%   primary:    axis label, passed to xlabel
%   unit:       unit label -- MATLAB figures out exponent and appends unit
function xlabels(primary, unit)
    arguments
        primary (1,1) string;
        unit (1,1) string;
    end
    import lplot.utl.get_ax_text;

    xaxis = gca().XAxis;
    unit = get_ax_text(xaxis, unit);
    
    xlabel(primary);
    xsecondarylabel(unit);
end
