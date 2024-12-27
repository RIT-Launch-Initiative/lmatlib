% NOTE: Requires Mapping Toolbox, <nwpdata>, and <xarray>
classdef atmosphere
    properties (GetAccess = public, SetAccess = private)
        data xarray;
        raster;
    end

    properties (Access = public)
        intrinsic_lat (1, 1) string;
        intrinsic_lon (1, 1) string;
        interpolant (1, 1) griddedInterpolant;
    end
    
    methods
        function obj = atmosphere(paths, params)
            arguments
                paths (:, 1) {mustBeFile};
                params.lat (1,1) double = 0;
                params.lon (1,1) double = 0;
                params.side (1,1) double = Inf;
            end
            
            data_array = cell(1, length(paths));
            for i = 1:length(paths)
                [data_array{i}, raster] = atmosphere.read_baro(paths, ...
                    elements = ["HGT", "UGRD", "VGRD", "TMP", "TKE"], ...
                    lat = params.lat, lon = params.lon, side = params.side);
            end

            obj.raster = raster;
            raster_type = obj.raster.CoordinateSystemType;
            if raster_type == "planar"
                obj.intrinsic_lat = "y";
                obj.intrinsic_lon = "x";
            elseif raster_type == "geographic"
                obj.intrinsic_lat = "lat";
                obj.intrinsic_lon = "lon";
            end

            data = cat("time", data_array{:});
            data = permute(data, [obj.intrinsic_lat, obj.intrinsic_lon, "time", "layer", "element"]);
            data.time = data.time - data.time(1);
            data = sort(data, [obj.intrinsic_lat, obj.intrinsic_lon, "time"], "ascend");
            obj.data = data;

            if isscalar(data.time)
                obj.interpolant = griddedInterpolant(data.coordinates(1:2), double(data), ...
                    "linear", "none");
            else
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

        function out = sample(obj, time, lat, lon)
            [int_lat, int_lon] = obj.get_intrinsic(lat, lon);
            if length(obj.interpolant.GridVectors) == 2
                out = obj.interpolant(int_lat, int_lon);
            elseif length(obj.interpolant.GridVectors) == 3
                out = obj.interpolant(int_lat, int_lon, seconds(time));
            end
            out = xarray(squeeze(out), layer = obj.data.layer, element = obj.data.element);
        end
        
        function out = tabulate(obj, time, lat, lon)
            samp = obj.sample(time, lat, lon);
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

        function [data, raster, metadata] = read_baro(path, params)
            % Read barometric level data from GRIB2 file into xarray
            % [data, raster, metadata] = atmosphere.read_baro(path, Name = Value)
            %   INPUTS
            %   path        Path to GRIB2 file
            %   
            %   Required name-value arguments
            %   elements    List of data elements (TMP, HGT, ...)
            %   
            %   Optional name-value arguments are used to crop the read data to
            %   a square to save memory. For planar projected maps, lat and lon
            %   are projected onto the map to get y/x limits. For geographic
            %   maps, side() is converted to angle limits using the Earth's
            %   radius, and the mean of the latitude limits is used to find the
            %   angle limits in longitude
            %   lat         Latitude enter of crop
            %   lon         Longitude of center of crop
            %   side        Side length of square used to crop
            %   
            %   OUTPUTS
            %   data        xarray with dimensions (y/x or lat/lon, layer, element, time)
            %               layer is pressure in Pa, sorted in descending order
            %   raster      Location referencing object
            %   metadata    Table containing information about layers and elements
            arguments
                path (1,1) string {mustBeFile};
                params.elements (1,:) string;
                params.lat (1, 1) double = 0;
                params.lon (1, 1) double = 0;
                params.side (1, 1) double = Inf;
            end

            isbl_regex = "(\d+)-ISBL";
            [data, raster, metadata] = nwpdata.read(path, ...
                elements = params.elements, layers = isbl_regex, ...
                lat = params.lat, lon = params.lon, side = params.side);
            baro_levels = string(regexp(string(data.layer), isbl_regex, "tokens"));
            data.layer = str2double(baro_levels);
            data = sort(data, "layer", "descend");
        end
    end

    methods (Static, Access = private)
    end
end

