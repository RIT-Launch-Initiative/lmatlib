% NOTE: Requires Mapping Toolbox, <nwpdata>, and <xarray>
classdef atmosphere
    properties (GetAccess = public, SetAccess = protected)
        data % xarray;
        raster % {mustBeA(raster, ["map.rasterref.MapPostingsReference", ...
            % "map.rasterref.GeographicPostingsReference"])};
        epoch (1,1) {mustBeA(epoch, ["datetime", "duration"])} = seconds(0);
    end

    properties (Access = protected)
        intrinsic_lat (1, 1) string;
        intrinsic_lon (1, 1) string
        sample3_template (:, :) %xarray;
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
            isbl_regex = "(\d+)-ISBL";

            for i = 1:length(paths)
                [data, raster] = nwpdata.read(paths(i), ...
                    fields = params.fields, layers = isbl_regex, ...
                    lats = params.lats, lons = params.lons);
                baro_levels = string(regexp(string(data.layer), isbl_regex, "tokens"));
                data.layer = str2double(baro_levels);
                data = data.rename("layer", "pressure");
                data_array{i} = sort(data, "pressure", "descend");
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
            data = cat("time", data_array{:});
            data = permute(data, [obj.intrinsic_lat, obj.intrinsic_lon, "time", "pressure", "field"]);
            data = sort(data, [obj.intrinsic_lat, obj.intrinsic_lon, "time"], "ascend");
            obj.epoch = data.time(1);
            data.time = seconds(data.time - obj.epoch);

            if isscalar(data.time)
                obj.interpolant = griddedInterpolant(data.coordinates(1:2), double(data), ...
                    "linear", "nearest");
            else
                obj.interpolant = griddedInterpolant(data.coordinates(1:3), double(data), ...
                    "linear", "nearest");
            end
            obj.data = data;
            obj.sample3_template = data.index(...
                obj.intrinsic_lat, 1, obj.intrinsic_lon, 1, "time", 1).squeeze;
        end

        function [int_lat, int_lon, int_time] = get_intrinsic(obj, lat, lon, time)
            % [int_lat, int_lon, int_time] = get_intrinsic(atmos, lat, lon, time)
            % Convert geographic-datetime coordinates to intrinsic (interpolant) coordinates
            % Time is converted to seconds since epoch
            % Lat/lon are returned unmodified 
            raster_type = obj.raster.CoordinateSystemType;
            if raster_type == "planar"
                [int_lon, int_lat] = projfwd(obj.raster.ProjectedCRS, lat, lon);
            elseif raster_type == "geographic"
                int_lat = lat;
                int_lon = lon;
            end
            int_time = seconds(time - obj.epoch);
        end

        function result = sample3(obj, time, lat, lon)
            arguments
                obj atmosphere;
                time (1,1) {mustBeA(time, ["datetime", "duration"])};
                lat (1,1) double;
                lon (1,1) double;
            end
            if length(obj.interpolant.GridVectors) == 2
                [int_lat, int_lon, ~] = obj.get_intrinsic(lat, lon, obj.epoch);
                out = obj.interpolant(int_lat, int_lon);
            elseif length(obj.interpolant.GridVectors) == 3
                [int_lat, int_lon, int_time] = obj.get_intrinsic(lat, lon, time);
                out = obj.interpolant(int_lat, int_lon, int_time);
            end

            result = obj.sample3_template;
            result{:, :} = squeeze(out);
        end

        function varargout = sample4(obj, time, lat, lon, height, fields)
            arguments
                obj atmosphere;
                time (1,1) {mustBeA(time, ["datetime", "duration"])};
                lat (1,1) double;
                lon (1,1) double;
                height (1,1) double;
                fields (1, :) string;
            end

            if nargout ~= length(fields)
                warning("%d outputs requested, %d used", length(fields), nargout);
            end

            samp = sample3(obj, time, lat, lon);

            interp_data = [samp.pressure, samp.double];
            result = interp1(samp.pick(field = "HGT").squeeze.double, interp_data, height);
            result = xarray(result, field = ["PRES"; samp.field]);

            present = ismember(fields, result.field);
            if any(~present)
                error("Data fields %s not found", fields(present));
            end
            result = result.pick(field = fields).double;

            varargout = cell(1, length(fields));
            for i = 1:length(fields)
                varargout{i} = result(i);
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

