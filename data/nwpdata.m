% Create and process NOAA's Numerical Weather Prediction data products
% Requires Mapping Toolbox for READGEORASTER
% 
% Desired functionality
% - Read data geographically limited by (lat, lon)
% - 

classdef (Sealed) nwpdata
    properties (Constant)
        sfc_match = "0-SFC";
        isbl_match = "(\d+)-ISBL";
    end

    methods (Static)
        function [tab] = read_tables(info, elements, lat, lon)
            tab = table;
            for element = elements
                [data, levels, bands] = nwpdata.read_baro(info, element, lat, lon);          
                [levels, order] = sort(levels, "descend");
                tab_element = table(levels(:), data(order), ...
                    VariableNames = ["ISBL", element]);

                metadata = info.Metadata(bands(1), :);
                tab_element.Properties.VariableUnits = ["[Pa]", metadata.Unit];
                tab_element.Properties.VariableDescriptions = ["Isobaric surface", ...
                    metadata.Comment];

                if isempty(tab)
                    tab = tab_element;
                else
                    tab = innerjoin(tab, tab_element, Keys = "ISBL");
                end
            end
            
            % Assign table metadata
            [~, name, ~] = fileparts(info.Filename);
            grid_name = regexp(name, "([A-Za-z]+_\d+?)_", "tokens");
            if isempty(grid_name)
                warning("Could not identify grid name");
                name = "unknown";
            else
                name = grid_name{1};
            end

            tab.Properties.UserData = struct(lat = lat, lon = lon, ...
                time = metadata.ReferenceTime(1), ...
                grid = name);
        end

        function [data, bands] = read_surf(info, element, lat, lon)
            sfc_level = "0-SFC";
            metadata = info.Metadata;

            bands = nwpdata.find_bands(metadata, element, sfc_level);
            if isempty(bands)
                error("No bands found")
            end

            fprintf("Reading %d bands from raster %s...\n", length(bands), info.Filename);
            tic; 
            [data, raster] = readgeoraster(info.Filename, Bands = bands); 
            toc;
        end

        function [data, levels, bands] = read_baro(info, element, lat, lon)
            baro_level = "(\d+)-ISBL";
            
            metadata = info.Metadata;
            bands = nwpdata.find_bands(metadata, element, baro_level);
            if isempty(bands)
                error("No bands found")
            end

            metadata = metadata(bands, :);
            levels_str = string(regexp(metadata.ShortName, baro_level, "tokens"));
            levels = str2double(levels_str);


            fprintf("Reading %d bands from raster %s...\n", length(bands), info.Filename);
            tic; 
            [data, raster] = readgeoraster(info.Filename, Bands = bands); 
            toc;

            data = nwpdata.get_point(data, raster, lat, lon);
            data = squeeze(data);
        end

        function [data_point, lat, lon] = get_point(data, raster, lat, lon)
            [x_c, y_c] = projfwd(raster.ProjectedCRS, lat, lon);
            x_lims = x_c + raster.SampleSpacingInWorldX/2*[-1 1];
            y_lims = y_c + raster.SampleSpacingInWorldY/2*[-1 1];
            [data_point, ~] = mapcrop(data, raster, x_lims, y_lims);
            data_point = squeeze(data_point(2, 2, :));
        end

        function bands = find_bands(metadata, element, level)
            element_idx = ~cellfun(@isempty, regexp(metadata.Element, element));
            level_idx = ~cellfun(@isempty, regexp(metadata.ShortName, level));
            bands = find(element_idx & level_idx);
        end
    end
end
