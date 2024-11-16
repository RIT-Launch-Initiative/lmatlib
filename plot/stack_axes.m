%% Stacks axes in specified TiledLayout object
% stack_axes(layout)
% - Narrows spacing
% - Gets rid of all ticks and tick labels except the last
% - Links all axes to the last
% 
% % Example
% lay = tiledlayout(3,1);
% nexttile; plot(table.Time, table.Data1);
% nexttile; plot(table.Time, table.Data2);
% nexttile; plot(table.Time, table.Data3);
% stack_axes(lay);

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
