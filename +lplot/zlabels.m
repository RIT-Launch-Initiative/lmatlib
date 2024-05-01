% Add primary and secondary Z-axis labels to a plot
% zlabels(primary, unit)
%   primary:    axis label, passed to zlabel
%   unit:       unit label -- MATLAB figures out exponent and appends unit
function zlabels(primary, unit)
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
