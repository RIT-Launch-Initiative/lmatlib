function text = get_ax_text(ax, unit)
    arguments
        ax handle;
        unit (1,1) string;
    end

    exponent = ax.Exponent;
    if (exponent == 0)
        text = sprintf("%s", unit);
    else
        text = sprintf("\\times10^{%d} %s", exponent, unit);
    end
end
