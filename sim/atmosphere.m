% NOTE: Requires Mapping Toolbox, <nwpdata>, and <xarray>
classdef atmosphere < handle
    properties (GetAccess = public, SetAccess = protected)
        ref; % nwpdata
        data; % grid of data
        pressures (1,:) double;
        sampler function_handle;
        
        epoch (1,1) datetime = missing;
    end

    properties (Access = protected)
        method (1, 1) string;
        extrap (1, 1) string;
        istimed (1, 1) logical = false;
    end

    properties (Constant, Access = protected)
        atmos_fields = ["HGT", "UGRD", "VGRD", "TMP"];
        sfc_layer = "0-SFC"
        htgl_layers = ["10-HTGL", "80-HTGL"];
        prs_regex = "(\d+)(?=-ISBL)";
    end
    
    methods
        function obj = atmosphere(ref, params)
            arguments
                ref (1,1) nwpdata;
                params.lats (1, 2) double = [-Inf Inf];
                params.lons (1, 2) double = [-Inf Inf];
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
            % to cull to maximum/minimum pressure levels

            data_args = {"fields", atmosphere.atmos_fields, ...
                "lats", params.lats, "lons", params.lons};

            warning("off", "nwpdata:missingFields");
            data = obj.ref.read(data_args{:}, ...
                layers = [atmosphere.htgl_layers, atmosphere.prs_regex]);
            warning("on", "nwpdata:missingFields");

            sfc = obj.ref.read(data_args{:}, layers = "0-SFC", fields = "HGT");
            hgts = xarray([10 80], layer = atmosphere.htgl_layers, field = "HGT");
            sfc = sfc.squeeze;
            data.pick(layer = atmosphere.htgl_layers, field = "HGT") = sfc + hgts;

            matches = regexp(data.layer, atmosphere.prs_regex, "match");
            obj.pressures = NaN(length(data.layer), 1);
            for i_layer = 1:length(data.layer)
                if ~isempty(matches{i_layer})
                    obj.pressures(i_layer) = str2double(matches{i_layer});
                end
            end

            % Convert text layers to numeric pressures
            % pressures = regexp(data.layer, atmosphere.prs_regex, "tokens");
            % pressures = str2double(string(pressures));
            % data.layer = pressures;
            % data = data.rename("layer", "pressure");
            % data = sort(data, "pressure", "descend");

            % Convert time
            obj.epoch = data.time(1);
            data.time = seconds(data.time - obj.epoch);

            obj.method = params.method;
            obj.extrap = params.extrap;
            if isscalar(data.time)
                obj.istimed = false;
                axes = [vertical horizontal];
            else
                obj.istimed = true;
                axes = [vertical horizontal "time"];
            end
            obj.data = data;
            obj.sampler = data.interpolant(axes, obj.method, obj.extrap);
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

            aircolumn = obj.sampler(sample);
            aircolumn = aircolumn.squeeze.permute(["layer", "field"]);
            % swap height and pressure
            aircolumn.layer = aircolumn.pick(field = "HGT").double;
            aircolumn = aircolumn.rename("layer", "height");
            aircolumn.pick{"field", "HGT"} = obj.pressures;
            aircolumn.field(aircolumn.field == "HGT") = "PRES";
            aircolumn = sort(aircolumn, "height", "ascend");

            if isfinite(params.height)
                result = interp1(aircolumn.height, double(aircolumn), params.height);
                result = xarray(result, field = aircolumn.field, height = params.height);
            else
                result = aircolumn;
            end

            if params.windout == "native"
                return;
            end

            result = result.permute(["height", "field"]);
            u_idx = find("UGRD" == result.field);
            v_idx = find("VGRD" == result.field);

            if isa(obj.ref.crs, "projcrs")  
                % get the wind data and make it the right shape
                nout = size(result, "height");
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
                h = result.height;
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

