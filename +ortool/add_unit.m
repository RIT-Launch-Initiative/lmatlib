function add_unit(axis, unit)
    arguments
        axis (1,1) string {is_axis};
        unit (1,1) string;
    end

    ax = gca();
    ax = ax.(strcat(upper(axis), "Axis"));

    exponent = ax.Exponent;
    if (exponent == 0)
        text = sprintf("""%s""", unit);
    else
        text = sprintf("""\\times10^{%d} %s""", exponent, unit);
    end
    
    f_call = sprintf("%ssecondarylabel(%s)", axis, text);
    eval(f_call);
end

function is_axis(axis)
    if ~ismember(axis, ["x", "y", "z"])
        error("%s is not an axis", axis);
    end
end
