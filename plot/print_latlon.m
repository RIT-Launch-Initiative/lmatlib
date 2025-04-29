function [str] = print_latlon(lat, lon, digits)
    % Print signed decimal latitude and longitude with N/S and E/W to required number of digits
    % str = print_latlon(lat, lon[, digits])
    % Inputs
    %   lat     (double)            [deg] Latitude
    %   lon     (double)            [deg] Longitude
    %   digits  ({4} | double)      Number of digits to print (1 = 10-km accuracy)
    arguments
        lat (1, 1) double;
        lon (1, 1) double;
        digits (1, 1) double {mustBeInteger} = 4;
    end
    if lat >= 0
        lat_letter = "N";
    else
        lat_letter = "S";
    end

    if lon >= 0
        lon_letter = "E";
    else
        lon_letter = "W";
    end
    deg = char(176);
    fmt = sprintf("%%#.%df", digits);
    str = sprintf("%s%s %s%s", num2str(abs(lat), fmt), deg + lat_letter, ...
        num2str(abs(lon), fmt), deg + lon_letter);
end

