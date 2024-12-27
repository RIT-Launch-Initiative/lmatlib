clear;

%% Grid-projected weather data
set(groot, "defaultAxesNextPlot", "add");
baro_inventory = atmosphere.baroproducts;

date = datetime(2024, 12, 1);
date.TimeZone = "UTC";
forecast = hours(0);

urrg = launchsites("urrg");
nam_file = nwpdata.download("data", "nam", baro_inventory.nam_12km, date, forecast);

[sfc, raster, meta] = nwpdata.read(nam_file, elements = ["HGT", "TMP"], layers = "0-SFC", ...
    lat = urrg.lat, lon = urrg.lon, side = 200e3); 

atmos = atmosphere(nam_file, lat = urrg.lat, lon = urrg.lon, side = 200e3);

urrg_data = atmos.sample(hours(0), urrg.lat, urrg.lon);

figure(name = "Air column at selected point");
tiledlayout(1,4)

H = double(urrg_data.subset(element = "HGT"));

[isa_T, ~, isa_P, isa_rho] = atmosisa(H);

% Compare temperatures
nexttile;
plot(isa_T, H/1000, DisplayName = "ISA");
plot(273.15 + double(urrg_data.subset(element = "TMP")), H/1000, DisplayName = "Local");
yline(urrg.alt/1000, "--k", HandleVisibility = "off");
legend;
ylabel("Geopotential height [km]");
xlabel("Temperature [K]");

% Compare pressures
nexttile;
plot(isa_P/1000, H/1000, DisplayName = "ISA");
plot(urrg_data.layer/1000, H/1000, DisplayName = "Local");
plot(urrg_data.layer/1000, pressalt("km", urrg_data.layer, "Pa"), DisplayName = "Indicated");
yline(urrg.alt/1000, "--k", HandleVisibility = "off");
legend;
xlabel("Pressure [kPa]");

% Calculate and plot density
nexttile;
urrg_rho = (1/287.05) * urrg_data.layer' ./ (273.15 + double(urrg_data.subset(element = "TMP")));
plot(isa_rho, H/1000, DisplayName = "ISA");
plot(urrg_rho, H/1000, DisplayName = "Local");
yline(urrg.alt/1000, "--k", HandleVisibility = "off");
legend;
xlabel("Density [kg/m^3]");

% Wind components
nexttile;
plot(urrg_data.subset(element = "UGRD"), H / 1000, DisplayName = "Eastward");
plot(urrg_data.subset(element = "VGRD"), H / 1000, DisplayName = "Northward");
comps = double(urrg_data.subset(element = ["UGRD", "VGRD"]));
magn = vecnorm(comps, 2, 2);
plot(magn, H / 1000, DisplayName = "Magnitude");
legend;
xlabel("Wind speed [m/s]");

sgtitle(sprintf("Atmospheric profile at %s\n%s", ...
    print_latlon(urrg.lat, urrg.lon, 4), string(date, "yyyy MMMM dd HH:mm z")));

% Wind infomation
figure(name = "Terrain and wind grid");
[x, y] = ndgrid(sfc.x, sfc.y);
z_0 = double(permute(squeeze(sfc.subset(element = "HGT")), ["x", "y"]));

surf(x, y, z_0, ...
    EdgeColor = "none");

windinfo = squeeze(atmos.data.subset(element = ["UGRD", "VGRD", "HGT"])); % get wind data, ignore time
windinfo = windinfo.subset(layer = (windinfo.layer > 70e3) & (mod(windinfo.layer, 5e3) == 0)); % crop and decimate layer
windinfo = windinfo.subset(x = 1:2:length(windinfo.x), y = 1:2:length(windinfo.y)); % decimate position
windinfo = permute(windinfo, ["x", "y", "layer", "element"]);

[x, y, ~] = ndgrid(windinfo.x, windinfo.y, windinfo.layer);
quiver3(x, y, double(windinfo.subset(element = "HGT")), ...
    double(windinfo.subset(element = "UGRD")), ...
    double(windinfo.subset(element = "VGRD")), ...
    zeros(size(x)));

colorbartext("Ground level [m]");
view(3);
% axis("tight", "equal");

%% Geographically sampled weather data
conus = {"lat", 40, "lon", -99, "side", 5000e3};
gfs_file = nwpdata.download("data", "gfs", baro_inventory.gfs_p50, date, forecast);
[sfc, raster, metadata] = nwpdata.read(gfs_file, conus{:}, ...
    elements = ["HGT", "TMP"], layers = "0-SFC");

atmos = atmosphere(gfs_file, conus{:});

figure(name = "CONUS elevation model");


sfc = sort(sfc, ["lat", "lon"], "ascend");
z = double(sfc.subset(element = "HGT"));
exaggeration = 100;
[lat, lon] = ndgrid(sfc.lat, sfc.lon);
[latlim, lonlim] = geoquadline(lat, lon);
worldmap(double(sfc), raster);
surfm(latlim, lonlim, z, exaggeration*z);
colorbartext("Ground level [m]");

windinfo = squeeze(atmos.data.subset(element = ["UGRD", "VGRD", "HGT"])); % get wind data, ignore time
windinfo = windinfo.subset(layer = (windinfo.layer >= 40e3) & (mod(windinfo.layer, 20e3) < 1)); % crop and decimate layer
windinfo = windinfo.subset(lat = 1:3:length(windinfo.lat), lon = 1:3:length(windinfo.lon)); % decimate position
windinfo = permute(windinfo, ["lat", "lon", "layer", "element"]);

[lat, lon] = ndgrid(windinfo.lat, windinfo.lon);

deg_per_ms = 0.05; % quiverm automatic scaling is worthless, use this to scale 1 m/s wind speed to some degrees lat/lon
for lay = 1:length(windinfo.layer)
    quiver3m(lat, lon, exaggeration * double(windinfo.subset(element = "HGT", layer = lay)), ...
        deg_per_ms * double(windinfo.subset(element = "VGRD", layer = lay)), ...
        deg_per_ms * double(windinfo.subset(element = "UGRD", layer = lay)), ...
        zeros(size(lat)), "y", 0);
end


view(3);

