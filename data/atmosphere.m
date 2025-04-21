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
    messages = ref.search(ref(1).inventory, filters);
    layers = unique(ref(1).inventory.layer(messages));

    % convert <xxx> mb to numbers and compare
    pressure_levels = str2double(extract(layers, digitsPattern));
    layers_in_range = layers(pressure_levels > opts.minpres);
    
    data = ref.read_point(lat, lon, layer = layers_in_range, field = opts.fields);

    % convert time axis to numeric (seconds since first time) so interp() works
    epoch = data.time(1);
    data.time = seconds(data.time - epoch);
    data = interp(data, time = seconds(time - epoch));
    data.layer = str2double(extract(data.layer, digitsPattern));
    data = data.align(["layer", "field"]).sort("layer", "descend"); % layer should be column so it fits in a table

    atmos = array2table([Pa_per_mbar * data.layer, double(data)], ...
        VariableNames = ["PRES", data.field']);
end
