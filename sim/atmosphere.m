% NOTE: Requires Mapping Toolbox, <nwpdata>, and <xarray>
classdef atmosphere < handle
    properties (GetAccess = public, SetAccess = protected)
        ref % nwpdata;
        prs_data % xarray;
        prs_sampler function_handle;
        epoch (1,1) datetime = missing;
    end

    properties (Access = protected)
        method (1, 1) string;
        extrap (1, 1) string;
        % vertical (1, 1) string ...
        %     {mustBeMember(vertical, ["y", "lat"])} = "lat";
        % horizontal (1, 1) string ...
        %     {mustBeMember(horizontal, ["x", "lon"])} = "lon";
        istimed (1, 1) logical = false;
    end

    properties (Constant, Access = protected)
        atmos_fields = ["HGT", "UGRD", "VGRD", "TMP"];
        prs_regex = "(\d+)-ISBL";
    end
    
    methods
        function obj = atmosphere(ref, params)
            arguments
                ref (1,1) nwpdata;
                params.lats (1, 2) double = [-Inf Inf];
                params.lons (1, 2) double = [-Inf Inf];
                % params.prs (1, 2) double = 
                params.method (1, 1) string = "linear";
                params.extrap (1, 1) string = "none";
            end

            if ~isdownloaded(ref)
                error("NWP data not downloaded");
            end
            obj.ref = ref;
            if isa(obj.ref.crs, "geocrs")
                vertical = "lat";
                horizontal = "lon";
            elseif isa(obj.ref.crs, "projcrs")
                vertical = "y";
                horizontal = "x";
            end

            % TODO restrict read data using raster metadata
            prs_data = obj.ref.read(fields = atmosphere.atmos_fields, ...
                layers = atmosphere.prs_regex, lats = params.lats, lons = params.lons);

            % Convert text layers to numeric pressures
            pressures = regexp(prs_data.layer, atmosphere.prs_regex, "tokens");
            pressures = str2double(string(pressures));
            prs_data.layer = pressures;
            prs_data = prs_data.rename("layer", "pressure");
            prs_data = sort(prs_data, "pressure", "descend");

            % Convert time
            obj.epoch = prs_data.time(1);
            prs_data.time = seconds(prs_data.time - obj.epoch);

            obj.method = params.method;
            obj.extrap = params.extrap;
            if isscalar(prs_data.time)
                obj.istimed = false;
                axes = [vertical horizontal];
            else
                obj.istimed = true;
                axes = [vertical horizontal "time"];
            end
            obj.prs_data = prs_data;
            obj.prs_sampler = prs_data.interpolant(axes, obj.method, obj.extrap);
        end

        function result = sample(obj, params)
            arguments
                obj atmosphere;
                params.lat (1,1) double {mustBeFinite};
                params.lon (1,1) double {mustBeFinite};
                params.time (1,1) datetime = NaT;
                params.windout (1,1) string ...
                    {mustBeMember(params.windout, ["native", "geographic"])} = "native";
                params.height (1,1) = NaN;
            end

            if isa(obj.ref.crs, "projcrs")
                [x, y] = projfwd(obj.ref.crs, params.lat, params.lon);
                sample = {y, x};
            else
                sample = {params.lat, params.lon};
            end

            if obj.istimed
                if ~isfinite(params.time)
                    error("Sample time required for this atmosphere model")
                end
                sample{3} = seconds(params.time - obj.epoch);
            end

            aircolumn = obj.prs_sampler(sample);
            aircolumn = aircolumn.squeeze.permute(["pressure", "field"]);

            if isfinite(params.height)
                data = [aircolumn.pressure, double(aircolumn)];
                data_height = aircolumn.pick(field = "HGT").double;
                data = interp1(data_height, data, params.height);
                % use already-allocated xarray - may be marginally more performant
                % result = aircolumn.index(pressure = 1);
                % result{:} = data(2:end);
                % result.pressure = data(1);
                result = xarray(data(2:end), ...
                    fields = aircolumn.field, pressure = data(1));
            else
                result = aircolumn;
            end

            if params.windout == "native"
                return;
            end

            u_idx = find("UGRD" == result.field);
            v_idx = find("VGRD" == result.field);

            if isa(obj.ref.crs, "projcrs")  
                % get the wind data and make it the right shape
                nout = size(result, "pressure");
                uvdata = result.index(field = [u_idx, v_idx]);

                assert(isequal(size(uvdata), [nout, 2]), ...
                    "Wrong shape input for array operation")

                winds = [x, y] + uvdata;
                [wlat, wlon] = projinv(obj.ref.crs, ...
                    winds.pick(field = "UGRD").double, ...
                    winds.pick(field = "VGRD").double);

                dlat = wlat - params.lat;
                dlon = wlon - params.lon;

                R = obj.ref.crs.GeographicCRS.Spheroid.SemimajorAxis;
                h = result.pick(field = "HGT").double;
                northward = deg2rad(dlat) .* (R + h);
                eastward = deg2rad(dlon) .* cosd(params.lat) .* (R + h);

                assert(isequal(size(northward), [nout 1]), ...
                    "Output N/E winds wrong shape for assignment")

                result{:, u_idx} = eastward;
                result{:, v_idx} = northward;
            end

            result.field(u_idx) = "EWND";
            result.field(v_idx) = "NWND";
        end
    end
end

