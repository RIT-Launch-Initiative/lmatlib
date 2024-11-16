%% Default axis grids to "on" or "off"
% set_grids("on"), set_grids("off")
function set_grids(mode)
    grids = compose("defaultAxes%sGrid", ["X", "Y", "Z"]);
    for grid = grids
            set(groot, grid, mode);
    end
end
