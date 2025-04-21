function alt = elevationquery(lat, lon)
    % Read altitude data from the OpenMeteo API
    % alt = elevationquery(lat, lon)
    % Inputs
    %   lat     (double)    [deg] latitude 
    %   lon     (double)    [deg] longitude
    % Outputs
    %   alt     (double)    [m] altitude MSL
    arguments
        lat (1,1) double;
        lon (1,1) double;
    end
    
    wopts = weboptions(ContentType = "json");
    url = sprintf("https://api.open-meteo.com/v1/elevation?latitude=%s&longitude=%s", ...
        compose("%f", lat).join(","), compose("%f", lon).join(","));
    ret = webread(url, wopts);
    alt = ret.elevation;
end
