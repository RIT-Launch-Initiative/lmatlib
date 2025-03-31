classdef atmosphere
    properties
        airdata xarray = xarray();
    end

    properties (Access = protected)
        epoch (1,1) datetime = NaT;
        airsample function_handle;
        crs {mustBeA(crs, ["geocrs", "projcrs"])} = geocrs(4326);
    end

    properties (Constant, Access = protected)
        % ncep_data_fields = ["HGT", "UGRD", "VGRD"];
        ncep_data_fields = ["HGT", "UGRD", "VGRD", "TMP"]; % what most users will want, but extra expense for us
        ncep_pressure_levels = regexpPattern("(^\d+) mb");
    end

    methods
        function atmos = atmosphere(air, crs)
            arguments
                air xarray;
                crs {mustBeA(crs, ["geocrs", "projcrs"])};
            end

            if isa(crs, "geocrs")
                mustHaveAxes(air, ["lat", "lon"]);
                vertical = "lat";
                horizontal = "lon";
            elseif isa(crs, "projcrs")
                mustHaveAxes(air, ["y", "x"]);
                vertical = "y";
                horizontal = "x";
            end

            mustHaveAxes(air, ["field", "layer", "time"]);

            air = air.pick(field = atmosphere.ncep_data_fields);

            if isempty(air)
                error("Air data must have fields U/VGRD, HGT, TMP")
            end
            
            if isscalar(air.time)
                axes = [vertical horizontal];         
            else
                atmos.epoch = air.time(1);
                air.time = atmos.datetime2sec(air.time);
                axes = [vertical horizontal "time"];
            end

            atmos.crs = crs;
            atmos.airdata = air;
            atmos.airsample = air.interpolant(axes, "linear", "nearest");
        end

        function aircol = aircolumn(atmos, lat, lon, time)
            arguments
                atmos (1,1) atmosphere;
                lat (1,1) double;
                lon (1,1) double;
                time (1,1) datetime = NaT;
            end

            [v, h] = atmos.latlon2vh(lat, lon);
            if isfinite(atmos.epoch)
                if ~isfinite(time)
                    error("Time required")
                end
                t = atmos.datetime2sec(time);
                aircol = atmos.airsample(v, h, t);
            else
                aircol = atmos.airsample(v, h);
            end
            
            aircol = aircol.squeeze.permute(["layer", "field"]);
            [~, order] = sort(aircol.pick{"field", "HGT"});
            aircol = aircol.index(layer = order);

            if isa(atmos.crs, "projcrs")
                uv = aircol.pick{"field", ["UGRD", "VGRD"]};
                assert(size(uv, 2) == 2, "Wind vectors must be N by 2");

                jacob = projjacob(atmos.crs, "geographic", lat, lon);
                jacob = jacob ./ vecnorm(jacob, 2, 1);
                winds = jacob * uv';
                aircol.pick{"field", ["UGRD", "VGRD"]} = winds';

                % points = [h, v]
                % points = [h, v] + uv;
                % [wp_lat, wp_lon] = projinv(atmos.crs, points(:, 1), points(:, 2));
                %
                % dlat = wp_lat - lat;
                % dlon = wp_lon - lon;
                %
                % R = atmos.crs.GeographicCRS.Spheroid.SemimajorAxis;
                % h = aircol.pick{"field", "HGT"};
                %
                % u = deg2rad(dlon) .* cosd(lat) .* (R+h);
                % v = deg2rad(dlat) .* (R+h);
                % aircol.pick{"field", ["UGRD", "VGRD"]} = [u, v];
            end
        end

    end

    methods (Static)
        function atmos = from_ncep(refs, params)
            arguments
                refs (1,:) ncep {mustBeNonempty};
                params.lats (1,2) double = [-Inf Inf];
                params.lons (1,2) double = [-Inf Inf];
                params.press (1,2) double = [-Inf Inf];
            end

            template = refs.inventory;
            layers = template.layer(matches(template.layer, ...
                atmosphere.ncep_pressure_levels));
            layers_press = str2double(extract(layers, digitsPattern));
            mask = min(params.press) < layers_press & layers_press < max(params.press);
            layers = unique(layers(mask));

            refs.download(field = atmosphere.ncep_data_fields, ...
                layer = layers);

            [airdata, crs] = refs.read(field = atmosphere.ncep_data_fields, ...
                layer = layers, lats = params.lats, lons = params.lons);

            atmos = atmosphere(airdata, crs);
        end
    end

    methods (Access = protected)
        function sec = datetime2sec(atmos, dt)
            arguments
                atmos atmosphere
                dt datetime
            end

            sec = seconds(dt - atmos.epoch);
        end

        function [v, h] = latlon2vh(atmos, lat, lon)
            if isa(atmos.crs, "geocrs")
                v = lat;
                h = lon;
            elseif isa(atmos.crs, "projcrs")
                [h, v] = projfwd(atmos.crs, lat, lon);
            end
        end
    end
end

