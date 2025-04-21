function atmos = atmosphere(model, product, lat, lon, time, opts)
    arguments
        model (1,1) string;
        product (1,1) string;
        lat (1,1) double;
        lon (1,1) double;
        time (1,1) datetime;
        opts.fields (1,:) string;
        opts.reftime (1,1) datetime;
        opts.minpres (1,1) double;
    end


end
