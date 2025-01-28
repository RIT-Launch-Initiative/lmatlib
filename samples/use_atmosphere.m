clear;
% This covers reading weather data with some settings. 
% samples/use_xarray.m sufficiently covers what one would do with the returned data.

data_folder = pfullfile("samples", "data");
% Model cycle
cycle = datetime(2024, 12, 1, TimeZone = "UTC") + hours(12);
urrg = launchsites("urrg");
% lat/lon limits for contiguous US
lat_limits = [24.5 49.5];
lon_limits = [-66 -125];

gfsref = nwpdata("gfs", "0.50 deg", cycle, hours(0:3:3));
gfsref.download(data_folder);

namref = nwpdata("nam", "12 km", cycle, hours(0:1:1));
namref.download(data_folder);

ref = namref;
atmos = atmosphere(ref, lats = lat_limits, lons = lon_limits);
urrg_air = atmos.sample(lat = urrg.lat, lon = urrg.lon, ...
    time = ref.cycle + ref.forecast(1));
[~, ~, pres, ~] = atmosisa(urrg_air.height);

figure;
plot(urrg_air.pick(field = "PRES"), urrg_air.height, DisplayName = ref.model);
plot(pres, urrg_air.height, DisplayName = "ISA");
xlabel("Pressure [Pa]");
ylabel("Height [gpm]");

figure;
urrg_air_N = atmos.sample(lat = urrg.lat, lon = urrg.lon, ...
    time = ref.cycle + ref.forecast(1), windout = "geographic");
plot(urrg_air.pick(field=["UGRD", "VGRD"]), urrg_air.height, "--");
plot(urrg_air_N.pick(field=["EWND", "NWND"]), urrg_air_N.height, "-");
xlabel("Wind component [m/s]");
ylabel("Height [gpm]");
