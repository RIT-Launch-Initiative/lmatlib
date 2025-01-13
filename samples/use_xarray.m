clear;

data_file = pfullfile("samples", "data", "gfs_wind.mat");
load(data_file);
disp(gfs_wind);
pressures = extract(gfs_wind.layer, digitsPattern); % The layer is formatted as <pressure>-ISBL, extract this number
pressures = str2double(pressures); % convert to double
gfs_wind.layer = pressures; % reassign text "layer" axis to use pressure
gfs_wind = sort(gfs_wind, "layer", "descend"); % sort by descending pressure (increasing height)

% subset of data, lat/lon range of roughly the continental US
conus = gfs_wind.range(lat = [25 50], lon = [-125 -70]); 
mid_alt = conus.pickt(layer = [80e3 10e3]); % pick data 80 +/- 10 kPa 
high_alt = conus.range(layer = [-Inf 70e3]); % pick data up to 70 kPa
mid_and_high_alt = cat("layer", mid_alt, high_alt); % glue data back together using dimension name

roc = mid_alt.pick(lat = 43, lon = -77.5);
h = roc.pick(element = "HGT");
u = roc.pick(element = "UGRD");
v = roc.pick(element = "VGRD");
mag = (u^2 + v^2)^0.5;

T_kelvin = roc.pick(element = "TMP") + 273.15;
rho = (1/287.05) * T_kelvin.layer / squeeze(T_kelvin); % re-arranged specific form of ideal gas law

figure;
hold on;
plot(mid_alt.layer, mid_alt.pick(element = "HGT", lat = 43, lon = -77)); % plot height against pressure level
plot(mid_and_high_alt.layer, mid_and_high_alt.pick(element = "HGT", lat = 43, lon = -77), ":");
xlabel("Pressure [Pa]");
ylabel("Height [gpm]");

figure;
hold on;
plot(u, h, DisplayName = "eastward");
plot(v, h, DisplayName = "northward");
plot(mag, h, DisplayName = "total");
xlabel("Wind speed [m/s]");
ylabel("Height [gpm]");
legend;

figure;
hold on;
plot(rho, h);
xlabel("Air density [kg/m^3]");
ylabel("Height [gpm]");
