function gh = plot3_openrocket(varargin)
    % Plot
    % line = plot3_openrocket([ax, ]data, ...)
    %   ax      ({gca} | Axes)  (Optional) Axes to target
    %   data    (timetable)     Output of openrocket.get_data() or import_openrocket_csv()
    %   ...     Additional arguments to plot ("-k", LineWidth=2, etc...)

    narginchk(1,Inf);

    if isgraphics(varargin{1}, "axes")
        ax = varargin{1};
        data = varargin{2};
        plot_args = varargin(3:end);
    elseif istimetable(varargin{1})
        ax = gca;
        view(ax, 3);
        data = varargin{1};
        plot_args = varargin(2:end);
    else
        error("First argument must be axis or timetable")
    end

    gh = plot3(data.("Position East of launch"), ...
        data.("Position North of launch"), ...
        data.("Altitude"), plot_args{:});
    
    if ax.XLabel.String == ""
        ax.XLabel.String = sprintf("East [%s]", ...
            data.Properties.VariableUnits{"Position East of launch"});
    end
    if ax.YLabel.String == ""
        ax.XLabel.String = sprintf("North [%s]", ...
            data.Properties.VariableUnits{"Position North of launch"});
    end
    if ax.ZLabel.String == ""
        ax.ZLabel.String = sprintf("Altitude [%s]", ...
            data.Properties.VariableUnits{"Altitude"});
    end
    if ax.DataAspectRatioMode == "auto"
        ax.DataAspectRatio = [1 1 1];
    end
end
