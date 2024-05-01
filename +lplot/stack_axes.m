function stack_axes(layout)
    layout.TileSpacing = "tight";
    layout.TileIndexing = "columnmajor";
    n = layout.GridSize(1);
    m = layout.GridSize(2);

    for col = 1:m
        last = nexttile(layout, n * col);

        for row = 1:n-1
            ax = nexttile(layout, n * (col-1) + row);

            % get rid of redundant texxt
            xticklabels(ax, {});
            xticklabels(ax, "manual");
            xlabel(ax, "");
            xsecondarylabel(ax, "");

            linkaxes([last ax], "x");
        end
    end
end
