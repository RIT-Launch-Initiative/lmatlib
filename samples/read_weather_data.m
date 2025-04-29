clear;

% valid time - time we care about
validtime = datetime(2024, 06, 21, 10, 20, 00, TimeZone = -hours(6));
% reference time - time the forecast is published
reftime = datetime(2024, 06, 21, 10, 20, 00, TimeZone = -hours(6)) - hours(12);

%% Demonstrate basic reading and plotting of one variable
figure(name = "Surface temperature");
layout = tiledlayout(2,1);
sgtitle("Surface temperature analysis nearest " + ...
    string(validtime, "dd-MMM-yyyy HH:mm z"));

% basic reading and plotting
disp("Models available:");
disp(ncep.list);

% create a referencing object for an HRRR analysis closest to the time we care about
hrrr_ref = ncep.analysis("hrrr", "wrfprsf", validtime);
disp("HRRR data inventory:");
disp(hrrr_ref.inventory);
% pick surface temperature -- from reading the inventory, the field is "TMP" and layer is "surface"
[hrrr_tempdata, ~] = hrrr_ref.read(layer = "surface", field = "TMP");
% output is an xarray with the identifying metadata (x/y/time/field/layer)

nexttile;
% 2-D xarrays can be plotted on a color-mapped image
imagesc(hrrr_tempdata, clabel = "Temperature [C]");
daspect([1 1 1]);
title("HRRR surface temperature at " + string(hrrr_tempdata.time));

% demonstrate geographically limiting data using lat and lon limits
nam_ref = ncep.analysis("nam", "awphys", validtime);
[nam_tempdata, ~] = nam_ref.read(layer = "surface", field = "TMP", ...
    lats = [40.5 49.5], lons = [-75 -92]);

nexttile;
imagesc(nam_tempdata, clabel = "Temperature [C]");
title("NAM 12-km surface temperature near Great Lakes at " + string(nam_tempdata.time));
daspect([1 1 1]);

sam = launchsites("spaceport-america");

%% Plotting wind forecasts
% read forecast from North American Mesoscale
% now that valid-time argument is a two-element vector, forecast() returns the
% tighest series of times that surrounds the specified bounds
nam_ref = ncep.forecast("nam", "awphys", validtime + hours([-1 1]), reftime);

% ncep.read_point reads a specific lat/lon point, taking care of coordinate
% projection and wind un-projection for you
nam_wind = nam_ref.read_point(sam.lat, sam.lon, ...
    layer = digitsPattern + " mb", field = ["UGRD", "VGRD", "HGT", "TMP"]);
nam_wind = sort_pressure_levels(nam_wind);

% paren-indexing into xarray returns another xarray
nam_wind = nam_wind.index(time = 1).range(layer = [400 1000]).align("layer"); 
% pick first time point, pressure between 

% read forecast from Global Forecast System
gfs_ref = ncep.forecast("gfs", "pgrb2.0p25", validtime + hours([-1 1]), reftime);
gfs_wind = gfs_ref.read_point(sam.lat, sam.lon, ...
    layer = digitsPattern + " mb", field = ["UGRD", "VGRD", "HGT", "TMP"]);
gfs_wind = sort_pressure_levels(gfs_wind);

gfs_wind = gfs_wind.index(time = 1).range(layer = [400 1000]).align("layer");

figure(name = "Forecast wind column for Spaceport America");
layout = tiledlayout(1,3);

nexttile; grid on; hold on;
view(3);

% this time, when we pick we use braces {} to return doubles
zr = zeros(size(nam_wind, "layer"), 1);
quiver3(zr, zr, nam_wind.pick{"field", "HGT"}, ...
    nam_wind.pick{"field", "UGRD"}, nam_wind.pick{"field", "VGRD"}, zr, "off", ...
    DisplayName = "NAM");

zr = zeros(size(gfs_wind, "layer"), 1);
quiver3(zr, zr, gfs_wind.pick{"field", "HGT"}, ...
    gfs_wind.pick{"field", "UGRD"}, gfs_wind.pick{"field", "VGRD"}, zr, "off", ...
    DisplayName = "GFS");

daspect([1 1 100]);
xlabel("Eastward wind [m/s]");
ylabel("Northward wind [m/s]");
zlabel("Altitude [m]");
title("Wind Vectors")
legend;

nexttile; grid on; hold on;
title("Pressure/Altitude");

plot(nam_wind.layer, nam_wind.pick{"field", "HGT"}, DisplayName = "NAM")
plot(gfs_wind.layer, gfs_wind.pick{"field", "HGT"}, DisplayName = "GFS")
xlabel("Pressure [mbar]");
ylabel("Geopotential altitude [m]");

nexttile; grid on; hold on;
title("Temperature/Altitude");

plot(nam_wind.pick{"field", "TMP"}, nam_wind.pick{"field", "HGT"}, DisplayName = "NAM");
plot(gfs_wind.pick{"field", "TMP"}, gfs_wind.pick{"field", "HGT"}, DisplayName = "GFS");
legend;
xlabel("Temperature [C]");
ylabel("Geopotential altitude [m]");


function data = sort_pressure_levels(data)
    data.layer = str2double(extract(data.layer, digitsPattern));
    data = sort(data, "layer", "descend");
end
