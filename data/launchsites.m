function info = launchsites(launchsite)
    siteinfo = dictionary;
    siteinfo("spaceport-america") = struct(lat = 32.94030518807581, lon = -106.92182007667907, alt = 1381);
    siteinfo("spaceport-midland") = struct(lat = 31.9485582, lon = -102.2086869, alt = 875);
    siteinfo("cedar-springs") = struct(lat = 37.7481336, lon = -113.2718213, alt = 1667);
    siteinfo("urrg") = struct(lat = 42.702986317065545, lon = -77.19120332086317, alt = 271);
    siteinfo("mars") = struct(lat = 42.799271225061474, lon = -77.84172220025518, alt = 170);
    
    info = siteinfo(launchsite);
end
