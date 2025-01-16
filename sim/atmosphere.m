% NOTE: Requires Mapping Toolbox, <nwpdata>, and <xarray>
classdef atmosphere
    properties (GetAccess = public, SetAccess = protected)
        data;
        raster;
    end

    properties (Access = protected)
        intrinsic_lat (1, 1) string;
        intrinsic_lon (1, 1) string;
        interpolant (1, 1) griddedInterpolant;
    end
    
    methods
        function obj = atmosphere(paths, params)
            arguments
                paths (:, 1) {mustBeFile};
                params.fields (1, :) string = ["HGT", "UGRD", "VGRD", "TMP"];
                params.lats (1, 2) double = [-Inf Inf];
                params.lons (1, 2) double = [-Inf Inf];
            end

            data_array = cell(1, length(paths));
            times = NaT(1, length(paths));
            isbl_regex = "(\d+)-ISBL";

            for i = 1:length(paths)
                [data, raster, meta] = nwpdata.read(paths(i), ...
                    fields = params.fields, layers = isbl_regex, ...
                    lats = params.lats, lons = params.lons);
                baro_levels = string(regexp(string(data.layer), isbl_regex, "tokens"));
                data.layer = str2double(baro_levels);
                data = data.rename("layer", "pressure");

                data_array{i} = sort(data, "pressure", "descend");
                times(i) = meta.ValidTime(1);
            end

            obj.raster = raster;
            raster_type = obj.raster.CoordinateSystemType;
            % input case: planar or geographic
            if raster_type == "planar"
                obj.intrinsic_lat = "y";
                obj.intrinsic_lon = "x";
            elseif raster_type == "geographic"
                obj.intrinsic_lat = "lat";
                obj.intrinsic_lon = "lon";
            end

            % input case: many time or 1 time
            if isscalar(paths)
                data = data_array{1};
                data = permute(data, [obj.intrinsic_lat, obj.intrinsic_lon, "pressure", "field"]);
                data = sort(data, [obj.intrinsic_lat, obj.intrinsic_lon], "ascend");
                obj.data = data;
                obj.interpolant = griddedInterpolant(data.coordinates(1:2), double(data), ...
                    "linear", "none");
            else
                data = cat("time", data_array{:});
                data.time = times;
                data = permute(data, [obj.intrinsic_lat, obj.intrinsic_lon, "time", "pressure", "field"]);
                data = sort(data, [obj.intrinsic_lat, obj.intrinsic_lon, "time"], "ascend");
                obj.data = data;
                obj.interpolant = griddedInterpolant(data.coordinates(1:3), double(data), ...
                    "linear", "none");
            end
        end

        function [int_lat, int_lon] = get_intrinsic(obj, lat, lon)
            raster_type = obj.raster.CoordinateSystemType;
            if raster_type == "planar"
                [int_lon, int_lat] = projfwd(obj.raster.ProjectedCRS, lat, lon);
            elseif raster_type == "geographic"
                int_lat = lat;
                int_lon = lon;
            end
        end

        function out = sample3(obj, time, lat, lon)
            [int_lat, int_lon] = obj.get_intrinsic(lat, lon);
            if length(obj.interpolant.GridVectors) == 2
                out = obj.interpolant(int_lat, int_lon);
            elseif length(obj.interpolant.GridVectors) == 3
                out = obj.interpolant(int_lat, int_lon, seconds(time));
            end
            out = xarray(squeeze(out), pressure = obj.data.pressure, field = obj.data.field);
        end

        function varargout = sample4(obj, time, lat, lon, height, fields)
            samp = sample3(obj, time, lat, lon);

            interp_data = [samp.layer(:), samp.double];
            result = interp1(samp.pick(field = "HGT").double, interp_data, height);
            result = xarray(result, height = height, field = ["PRES", samp.field]);

            present = ismember(fields, result.field);
            if any(~present)
                error("Data fields %s not found", fields(present));
            end
            result = result.pick(field = fields);

            if nargout == 1
                varargout{1} = result;
            else
                unpack = mat2cell(result.double, length(result.height), ones(1, length(result.field)));
                [varargout{1:nargout}] = unpack{:};
            end
        end
    end

    methods (Static)
        function info = baroproducts
            info.gfs_100 = "pgrb2.1p00";
            info.gfs_p50 = "pgrb2.0p50";
            info.gfs_p25 = "pgrb2.0p25";
            info.nam_12km = "awphys";
            info.nam_3km = "conusnest.hiresf";
            info.hrrr_3km = "wrfprsf";
        end
    end
end

