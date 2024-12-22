classdef (Sealed) atmosphere
    properties (GetAccess = public, SetAccess = private)
        weathergrid nwpdata;
        x_grid (:, :, :) double;
        y_grid (:, :, :) double;
        t_grid (:, :, :) double;
        pressures (1, :) double;
    end
    
    methods
        function obj = atmosphere(dataset)
            baro_layers = regexp(dataset.layers, "(\d+)-ISBL", "tokens");
            baro_layers = str2double(string(baro_layers));
            [obj.pressures, l_order] = sort(baro_layers, "descend");
            [~, x_order] = sort(dataset.proj_x, "ascend");
            [~, y_order] = sort(dataset.proj_y, "ascend");

            obj.weathergrid = dataset(:, x_order, y_order, l_order, ["HGT", "UGRD", "VGRD", "TMP", "TKE"]);
            [obj.t_grid, obj.x_grid, obj.y_grid] = ...
                ndgrid(seconds(obj.weathergrid.time), obj.weathergrid.proj_x, obj.weathergrid.proj_y);
        end

        function [x, y] = project(obj, lat, lon)
            [x, y] = projcrs(obj.weathergrid.projection, lat, lon);
        end

        function [tabl] = tabulate(obj, t, x, y)
            samp = interpn(obj.t_grid, obj.x_grid, obj.y_grid,  ...
                obj.weathergrid.data, seconds(t), x, y);
            samp = squeeze(samp);
            tabl = array2table([obj.pressures(:), samp], ...
                VariableNames = ["PRES", obj.weathergrid.quantities]);
        end

        function varargout = sample(obj, t, x, y, z, q)
            samp_txy = interpn(obj.t_grid, obj.x_grid, obj.y_grid,  ...
                obj.weathergrid.data, seconds(t), x, y);
            samp_txy = [obj.pressures', squeeze(samp_txy)];
            samp_names = ["PRES", obj.weathergrid.quantities];
            hgt_index = samp_names == "HGT";
            [~, ~, queried] = intersect(q, samp_names, "stable");
            samp_txyzq = interp1(samp_txy(:, hgt_index), samp_txy(:, queried), z);

            varargout = cell(1, nargout);
            if nargout == length(q)
                for i = 1:nargout
                    varargout{i} = samp_txyzq(i);
                end
            else
                varargout{1} = samp_txyzq;
            end
            
        end
    end
end

