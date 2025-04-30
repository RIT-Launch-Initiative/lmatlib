% Create atmospheric data table using a specific model, grid, and space-time
% atmos = atmosphere(model, product, lat, lon, time, [Name, Value])
% INPUTS
%   model   (string)    NOAA model to pull from - use >>ncep.list()
%   product (string)    Model output grid to pull from - use >>ncep.list(model)
%   lat     (double)    Signed latitude, decimal degrees North
%   lon     (double)    Signed lonitude, decimal degrees East
% NAME-VALUE ARGUMENTS
%   fields  ({["TMP", "HGT", "UGRD", "VGRD"]} | string)
%           Data fields to pull, defaults to the typical wanted ones
%   reftime ({NaT} | datetime)
%           Reference time for forecast (time at which the forecast is generated)
%           Should have time zone assigned for minimum ambiguity
%           Rounds to the nearest model output cycle
%   minpres ({-Inf} | double)
%           Minimum pressure to reduce the number of layers (GRIB messages) that need to be read
% OUTPUTS
%   atmos   (table)     Table with columns ["PRES", fields] describing the
%                       atmosphere profiles per pressure level, up to minpres
function atmos = atmosphere(model, product, lat, lon, time, opts)
    arguments
        model (1,1) string;
        product (1,1) string;
        lat (1,1) double;
        lon (1,1) double;
        time (1,1) datetime;
        opts.fields (1,:) string = ["TMP", "HGT", "UGRD", "VGRD"];
        opts.reftime (1,1) datetime = missing;
        opts.minpres (1,1) double = -Inf;
        opts.cache (1,1) string = missing;
    end
    
    if ~ismissing(opts.cache)
        mat = matfile(opts.cache, Writable = true);
        cache_key = keyHash({model, product, lat, lon, time, ...
            opts.fields, opts.reftime, opts.minpres});
        cache_name = key2name(cache_key);
        if ~isempty(whos(mat, cache_name))
            atmos = mat.(cache_name);
            return;
        end
    end

    % Magic values
    Pa_per_mbar = 100;
    pressure_level_regex = "^\d+ mb$";

    if ismissing(opts.reftime)
        ref = ncep.analysis(model, product, [time time]);
    else
        ref = ncep.forecast(model, product, [time time], opts.reftime);
    end

    % determine pressure levels to use based on minpres
    filters.field = opts.fields;
    filters.layer = regexpPattern(pressure_level_regex);

    % inventory has <xxx> mb layer codes, find them
    inv = ref(1).inventory;
    messages = ref.search(inv, filters);
    layers = unique(inv.layer(messages));
    notpresent = setdiff(opts.fields, inv.field(messages));
    if ~isempty(notpresent)
        warning("Fields %s not present in inventory", mat2str(notpresent));
    end

    % convert <xxx> mb to numbers and compare
    pressure_levels = str2double(extract(layers, digitsPattern));
    layers_in_range = layers(pressure_levels > opts.minpres);
    
    data = ref.read_point(lat, lon, layer = layers_in_range, field = opts.fields);

    % if the specified result is on-the-hour, size(data, "time") = 1 and the interpolation fails
    if size(data, "time") > 1
        % convert time axis to numeric (seconds since first time) so interp() works
        epoch = data.time(1);
        data.time = seconds(data.time - epoch);
        data = interp(data, time = seconds(time - epoch));
    end
    data.layer = str2double(extract(data.layer, digitsPattern));
    data = data.align(["layer", "field"]).sort("layer", "descend"); % layer should be column so it fits in a table

    atmos = array2table([Pa_per_mbar * data.layer, double(data)], ...
        VariableNames = ["PRES", data.field']);

    if ~ismissing(opts.cache)
        mat.(cache_name) = atmos;
    end
end

function name = key2name(key)
    name = sprintf("data_%ud", key);
end
