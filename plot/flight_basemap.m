function [im, rast] = flight_basemap(lat, lon, width)
    % Read a satellite basemap image and shift the returned raster to place
    % (lat, lon) at (0, 0) for flight data plotting. Requires Mapping Toolbox.
    %
    % [im, rast] = flightBasemapImage(lat, lon, width)
    %   lat     (double)    [deg N] Center point latitude
    %   lon     (double)    [deg E] Center point longitude
    %   width   (double)    [m] Map side length
    %
    %   im      (3-D double)            Color data
    %   rast    (MapPostingsReference)  Raster for plotting
    %
    % mapshow(im, rast) % plots map raster
    arguments
        lat (1,1) double;
        lon (1,1) double;
        width (1,1) double;
    end

    R = 6370e3;
    
    dlat = rad2deg(width / R) / 2;
    dlon = rad2deg(width / R / cosd(lat)) / 2;
    [im, rast] = readBasemapImage("satellite", lat + dlat * [-1 1], lon + dlon * [-1 1]);
    [x, y] = projfwd(rast.ProjectedCRS, lat, lon);
    rast.XWorldLimits = rast.XWorldLimits - x;
    rast.YWorldLimits = rast.YWorldLimits - y;
end
