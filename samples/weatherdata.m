clear;
set(groot, "defaultAxesNextPlot", "add");

%% Grid-projected weather data
baro_inventory = atmosphere.baroproducts;

date = datetime(2024, 12, 1);
date.TimeZone = "UTC";
forecast = hours(0);
urrg = launchsites("urrg");

profile on;
nam_file = nwpdata.download("data", "nam", baro_inventory.nam_12km, date, forecast);

[nam_surface, nam_raster, nam_meta] = nwpdata.read(nam_file, fields = ["HGT", "TMP"], ...
    layers = "0-SFC", ...
    lats = urrg.lat + [-0.5 0.5], lons = urrg.lon + [-0.5 0.5]); 

nam_atmos = atmosphere(nam_file, lats = urrg.lat + [-0.5 0.5], lons = urrg.lon + [-0.5 0.5]);

urrg_nam_data = nam_atmos.sample3(hours(0), urrg.lat, urrg.lon);

figure(name = "Air column at selected point");
tiledlayout(1,4)

H = urrg_nam_data.pick(field = "HGT");
profile viewer;

[isa_T, ~, isa_P, isa_rho] = atmosisa(double(H));

% Compare temperatures
nexttile;
plot(isa_T - 273.15, H/1000, DisplayName = "ISA");
plot(urrg_nam_data.pick(field = "TMP"), H / 1000, DisplayName = "Local");
yline(urrg.alt/1000, "--k", HandleVisibility = "off");
legend;
ylabel("Geopotential height [km]");
xlabel("Temperature [C]");

% Compare pressures
nexttile;
plot(isa_P/1000, H/1000, DisplayName = "ISA");
plot(urrg_nam_data.pressure/1000, H/1000, DisplayName = "Local");
plot(urrg_nam_data.pressure/1000, pressalt("km", urrg_nam_data.pressure, "Pa"), DisplayName = "Indicated");
yline(urrg.alt/1000, "--k", HandleVisibility = "off");
legend;
xlabel("Pressure [kPa]");

% Calculate and plot density
nexttile;
urrg_rho = (1/287.05) * urrg_nam_data.pressure / (273.15 + urrg_nam_data.pick(field = "TMP"));
plot(isa_rho, H/1000, DisplayName = "ISA");
plot(urrg_rho, H/1000, DisplayName = "Local");
yline(urrg.alt/1000, "--k", HandleVisibility = "off");
legend;
xlabel("Density [kg/m^3]");

% Wind components
nexttile;
plot(urrg_nam_data.pick(field = "UGRD"), H / 1000, DisplayName = "Eastward");
plot(urrg_nam_data.pick(field = "VGRD"), H / 1000, DisplayName = "Northward");
wind_magnitudes = (urrg_nam_data.pick(field = "UGRD").^2 + urrg_nam_data.pick(field = "UGRD").^2).^0.5;
plot(wind_magnitudes, H / 1000, DisplayName = "Magnitude");
legend;
xlabel("Wind speed [m/s]");

sgtitle(sprintf("Atmospheric profile at %s\n%s", ...
    print_latlon(urrg.lat, urrg.lon, 4), string(date, "yyyy MMMM dd HH:mm z")));

% Wind infomation
figure(name = "Terrain and wind grid");
[x, y] = ndgrid(nam_surface.x, nam_surface.y);
z_0 = nam_surface.pick(field = "HGT").permute(["x", "y"]).double;

surf(x, y, z_0, ...
    EdgeColor = "none");

% get wind (U, V, H) data
% take all pressures exceeding 70 kPa
% decimate x/y by 2 
% finally re-arrange
windinfo = nam_atmos.data;
windinfo = windinfo.pick(field = ["UGRD", "VGRD", "HGT"]) ...
    .range(pressure = [70e3 Inf]) ...
    .index(x = 1:2:size(windinfo, "x"), y = 1:2:size(windinfo, "y")) ...
    .permute(["x", "y", "pressure", "field"]);

[x, y, ~] = ndgrid(windinfo.x, windinfo.y, windinfo.pressure);
quiver3(x, y, windinfo.pick(field = "HGT").double, ...
    windinfo.pick(field = "UGRD").double, ...
    windinfo.pick(field = "VGRD").double, ...
    zeros(size(x)));

colorbartext("Ground level [m]");
view(3);
% axis("tight", "equal");

%% Geographically sampled weather data
conus = {"lats", [24.5 49.5], "lons", [-66 -125]};

gfs_file = nwpdata.download("data", "gfs", baro_inventory.gfs_p50, date, forecast);
[gfs_surface, gfs_raster, gfs_metasata] = nwpdata.read(gfs_file, conus{:}, ...
    fields = ["HGT", "TMP"], layers = "0-SFC");

gfs_atmos = atmosphere(gfs_file, conus{:});

figure(name = "CONUS winds at select layers");

gfs_surface = sort(gfs_surface, ["lat", "lon"], "ascend");
z = gfs_surface.pick(field = "HGT").double;
exaggeration = 100;

[lat, lon] = ndgrid(gfs_surface.lat, gfs_surface.lon);
[latlim, lonlim] = geoquadline(lat, lon);

worldmap(double(gfs_surface), gfs_raster); % initializes map-style axes using the geographic limits implied by gfs_raster
surfm(latlim, lonlim, z, exaggeration*z);
view(3);
colorbartext("Ground level [m]");

windinfo = gfs_atmos.data;
windinfo = windinfo.pick(field = ["UGRD", "VGRD", "HGT"]) ...
    .range(pressure = [40e3 Inf]) ...
    .index(lat = 1:2:size(windinfo, "lat"), lon = 1:2:size(windinfo, "lon")) ...
    .permute(["lat", "lon", "pressure", "field"]);
windinfo = windinfo.index(pressure = mod(windinfo.pressure, 20e3) == 0); % can't decimate until the .range on pressure change has applied

[lat, lon] = ndgrid(windinfo.lat, windinfo.lon);

% quiverm automatic scaling is worthless, use this to scale 1 m/s wind speed to some degrees lat/lon
deg_per_ms = 0.05; 
for lay = 1:length(windinfo.pressure)
    wind_at_layer = windinfo.index(pressure = lay);
    quiver3m(lat, lon, exaggeration * wind_at_layer.pick(field = "HGT").double, ...
        deg_per_ms * wind_at_layer.pick(field = "VGRD").double, ...
        deg_per_ms * wind_at_layer.pick(field = "UGRD").double, ...
        zeros(size(lat)), "y", 0);
end

urrg_gfs_data = gfs_atmos.sample3(hours(0), urrg.lat, urrg.lon);
profile viewer;

figure(name = "GFS vs NAM winds at selected location");
tiledlayout(1,2);

H_nam = urrg_nam_data.pick(field = "HGT");
H_gfs = urrg_gfs_data.pick(field = "HGT");

nexttile;
hold on;
plot(urrg_nam_data.pick(field = "UGRD"), H_nam / 1000, DisplayName = "NAM sample");
plot(urrg_gfs_data.pick(field = "UGRD"), H_gfs / 1000, DisplayName = "GFS sample");
legend;
xlabel("Eastward wind speed [m/s]");
ylabel("Geopotential height [km]");

nexttile;
hold on;
plot(urrg_nam_data.pick(field = "VGRD"), H_nam / 1000, DisplayName = "NAM sample");
plot(urrg_gfs_data.pick(field = "VGRD"), H_gfs / 1000, DisplayName = "GFS sample");
legend;
xlabel("Northward wind speed [m/s]");
