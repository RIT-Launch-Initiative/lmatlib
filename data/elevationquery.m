function alt = elevationquery(lat, lon, opts)
    % Read altitude data from the USGS Elevation Point Query Service
    % alt = elevationquery(lat, lon[, Name, Value])
    %   Inputs
    %   lat     (double)    latitude in degrees
    %   lon     (double)    longitude in degrees
    %   Name-Value arguments
    %   units   ({"m"} | "ft")  Output units
    %   Outputs
    %   alt     (double)    altitude
    % See documents for <a href="https://www.usgs.gov/faqs/how-accurate-are-elevations-generated-elevation-point-query-service-national-map">accuracy</a>
    arguments
        lat (1,1) double;
        lon (1,1) double;
        opts.units (1,1) string {mustBeMember(opts.units, ["ft", "m"])} = "m";
    end

    if opts.units == "m"
        unitstr = "Meters";
    elseif opts.units == "ft"
        unitstr = "Feet";
    end

    wopts = weboptions(ContentType = "json", Timeout = 10);
    basepath = "https://epqs.nationalmap.gov/v1/json";
    keys = sprintf("?x=%f&y=%f&units=%s&wkid=4326&includeDate=False", lon, lat, unitstr);
    ret = webread(basepath + keys, wopts);
    alt = ret.value;
    % Feet return numeric result, Meters return char result for some reason
    if ischar(alt) || isstring(alt)
        alt = str2double(alt);
    end
end

