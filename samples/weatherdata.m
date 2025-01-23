clear;
% This covers reading weather data with some settings. 
% samples/use_xarray.m sufficiently covers what one would do with the returned data.
%% Geographically sampled weather data
date = datetime(2024, 12, 1, TimeZone = "UTC");
data_folder = pfullfile("samples", "data");
conus = {"lats", [24.5 49.5], "lons", [-66 -125]};

gfsref = nwpdata("gfs", "0.50 deg", date, hours(0:3:6));
gfsref.download(data_folder);

sfctemp = gfsref.read(conus{:}, fields = "TMP", layers = "0-SFC");

%% Grid-projected weather data
namref = nwpdata("nam", "12 km", date, hours(0:1:5));
namref.download(data_folder);

isblheights = namref.read(fields = "HGT", layers = "(\d+)-ISBL");
