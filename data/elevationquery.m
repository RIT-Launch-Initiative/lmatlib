function alt = elevationquery(lat, lon)
    % alt = elevationquery(lat, lon)
    %   read altitude data from the USGS Elevation Point Query Service, in meters
    arguments
        lat (1,1) double;
        lon (1,1) double;
    end

    wopts = weboptions(ContentType = "json");
    basepath = "https://epqs.nationalmap.gov/v1/json";
    keys = sprintf("?x=%f&y=%f&units=Meters&wkid=4326&includeDate=False", lon, lat);
    ret = webread(basepath + keys, wopts);
    alt = str2double(ret.value);
end

