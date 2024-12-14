% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

classdef (Sealed) nwpdata
    properties
        data (:, :, :) double;
        raster map.rasterref.MapPostingsReference;
        levels (1, :) double;
        comment (1, 1) string;
    end

    properties (Constant)
        sfc_match = "0-SFC";
        isbl_match = "(\d+)-ISBL";
    end

    methods
        function obj = nwpdata(path, element, lat, lon, side)
            arguments
                path (1,1) string;
                element (1,1) string;
                lat (1,1) double;
                lon (1,1) double;
                side (1,1) double = 0;
            end

            tic;
            fprintf("Reading raster information...\n");
            info = georasterinfo(path);
            toc;

            meta = info.Metadata;
            raster = info.RasterReference;

            [x_c, y_c] = projfwd(raster.ProjectedCRS, lat, lon);
            if side == 0
                x_length = raster.SampleSpacingInWorldX;
                y_length = raster.SampleSpacingInWorldY;
            else
                x_length = side;
                y_length = side;
            end
            x_lims = x_c + x_length/2*[-1/2 1/2];
            y_lims = y_c + y_length/2*[-1/2 1/2];

            
            element = "^" + element + "$";
            level = "^" + nwpdata.isbl_match + "$";
            
            element_idx = ~cellfun(@isempty, regexp(meta.Element, element));
            level_idx = ~cellfun(@isempty, regexp(meta.ShortName, level));
            bands = find(element_idx & level_idx);
            if isempty(bands)
                error("No matching bands found");
            end

            tic;
            fprintf("Reading %d bands from %s.\n", length(bands), info.Filename)
            [data, rast] = readgeoraster(info.Filename, Bands = bands);
            toc;

            meta = meta(bands, :);
            [data, rast] = mapcrop(data, rast, x_lims, y_lims);
            levels_str = string(regexp(meta.ShortName, level, "tokens"));
            levels = str2double(levels_str);
            [levels, reorder] = sort(levels, "descend");
            
            obj.data = data(:, :, reorder);
            obj.raster = rast;
            obj.levels = levels;
        end
    end

    methods (Static)
        function [data, raster, levels, comment] = read_baro(info, element, x_lims, y_lims)
            level = nwpdata.isbl_match;
            bands = nwpdata.find_bands(element, level);
            
            if isempty(bands)
                error("No matching bands found");
            end

            tic;
            fprintf("Reading %d band(s) from %s.\n", length(bands), info.Filename)
            [data, rast] = readgeoraster(info.Filename, Bands = bands);
            toc;

            meta = meta(bands, :);
            comment = meta(1, "Comment");

            [data, rast] = mapcrop(data, rast, x_lims, y_lims);
            levels_str = string(regexp(meta.ShortName, "^" + level + "$", "tokens"));
            levels = str2double(levels_str);
            [levels, reorder] = sort(levels, "descend");

            data = data(:, :, reorder);
        end

        function [data, raster] = read_flat(element, level, x_lims, y_lims)
            bands = nwpdata.find_bands(element, level);
            if ~isscalar(bands)
                error("%d bands matching '%s' found, expected one", length(bands), level);
            end

            tic;
            fprintf("Reading %d band(s) from %s.\n", length(bands), info.Filename)
            [data, rast] = readgeoraster(info.Filename, Bands = bands);
            toc;

            meta = meta(bands, :);
            comment = meta.Comment;
            [data, rast] = mapcrop(data, rast, x_lims, y_lims);
        end

        function bands = find_bands(element, level)
            element = "^" + element + "$";
            level = "^" + level + "$";
            
            element_idx = ~cellfun(@isempty, regexp(meta.Element, element));
            level_idx = ~cellfun(@isempty, regexp(meta.ShortName, level));
            bands = find(element_idx & level_idx);
        end

        function dict = make_types
            dict("UGRD") = "wind-east";
            dict("VGRD") = "wind-north";
            dict("DZDT") = "wind-up";
            dict("ABSV") = "vorticity";
            dict("TMP") = "temperature";
            dict("HGT") = "height";
        end
    end
end
