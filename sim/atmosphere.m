% NOTE: Requires Mapping Toolbox, <nwpdata>, and <xarray>
classdef atmosphere < handle
    properties (GetAccess = public, SetAccess = protected)
        ref nwpdata;
        prs_data % xarray;
        prs_sampler function_handle;
        epoch (1,1) datetime = missing;
    end

    properties (Access = protected)
        method (1, 1) string;
        extrap (1, 1) string;
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

            result = obj.prs_sampler(sample).squeeze;
            if isfinite(params.height)
                data = [result.pressure, double(result)];
                data = interp1(result.height, data, params.height);
                result = xarray(data(2:end), ...
                    fields = result.fields, pressure = data(1));
            end
        end
    end
end

