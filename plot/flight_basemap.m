function [im, rast, attrib] = flight_basemap(lat, lon, width, opts)
    % Read a satellite basemap image and shift the returned raster to place
    % (lat, lon) at (0, 0) for flight data plotting. Requires Mapping Toolbox.
    %
    % [im, rast] = flight_basemap(lat, lon, width)
    % Inputs
    %   lat     (double)    [deg N] Center point latitude
    %   lon     (double)    [deg E] Center point longitude
    %   width   (double)    [m] Map side length
    % Name-Value Inputs
    %   cache   (string)    Name of .mat file storing basemaps
    % Outputs
    %   im      (3-D double)            Color data
    %   rast    (MapPostingsReference)  Raster for plotting
    %   attrib  (string)                Attribution
    %
    % mapshow(im, rast) % plots map raster
    arguments
        lat (1,1) double;
        lon (1,1) double;
        width (1,1) double;
        opts.cache (1,1) string = missing;
    end

    R = 6370e3;
    
    if ~ismissing(opts.cache)
        mat = matfile(opts.cache, Writable = true);
        cache_key = keyHash({lat, lon, width});
        cache_name = key2name(cache_key);

        if ~isempty(whos(mat, cache_name))
            data = mat.(cache_name);

            im = data.im;
            rast = data.rast;
            attrib = data.attrib;

            return;
        end
    end

    dlat = rad2deg(width / R) / 2;
    dlon = rad2deg(width / R / cosd(lat)) / 2;
    [im, rast, attrib] = readBasemapImage("satellite", lat + dlat * [-1 1], lon + dlon * [-1 1], 25);
    [x, y] = projfwd(rast.ProjectedCRS, lat, lon);
    rast.XWorldLimits = rast.XWorldLimits - x;
    rast.YWorldLimits = rast.YWorldLimits - y;

    if ~ismissing(opts.cache)
        data.im = im;
        data.rast = rast;
        data.attrib = attrib;

        mat.(cache_name) = data;
    end
end

function name = key2name(key)
    name = sprintf("data_%ud", key);
end
