function info = launchsites(launchsite)
    siteinfo = dictionary;

    siteinfo("spaceport-america") = struct(lat = 32.94012563817728, lon = -106.91165732518579, alt = 1381);
    siteinfo("spaceport-midland") = struct(lat = 31.049770791981764, lon = -103.54718796542103, alt = 875);
    siteinfo("urrg") = struct(lat = 42.702986317065545, lon = -77.19120332086317, alt = 271);
    siteinfo("mars") = struct(lat = 42.799271225061474, lon = -77.84172220025518, alt = 170);
    
    info = siteinfo(launchsite);
end
