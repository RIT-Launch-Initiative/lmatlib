clear;

load("data/gfswind.mat", "data", "axes", "coords");
atmosphere_data = xarray(data, axes, coords);

profile on;
%% Basic indexing and interpolation
% Get temperatures along the 43rd parallel at analysis-time
% use <pick> to select by matching values
% use <index> to select by axis index
%   index is like data(indices1,  indices2, ...), but you can specify a subset of the axes by name
temp_data = atmosphere_data.pick(lat = 43, field = "TMP").index(time = 1);
temp_at_pressure = temp_data.interpolant("pressure", "linear"); % create interpolation object

temp_near_surface = temp_data.pick(pressure = 100e3);
temp_aloft = temp_data.pick(pressure = 90e3);
temp_between = temp_at_pressure(96e3);

figure;
plot(temp_near_surface.lon, double(temp_near_surface), ...
    DisplayName = "1000 mbar layer"); % convert the data to a double for plotting
plot(temp_aloft.lon, double(temp_aloft), ...
    DisplayName = "900 mbar layer"); % convert the data to a double for plotting
plot(temp_aloft.lon, double(temp_between), ...
    DisplayName = "960 mbar interpolation"); % convert the data to a double for plotting
xlabel("Longitude [deg]");
ylabel("Temperature [C]");
legend;

% Exercise: modify the temp_data to only include approximately the range
% containing the continental U.S.

%% Interpolation
% Heights of the isobaric surfaces over approximately the contiguous U.S. at the 3-hour forecast
height_data = atmosphere_data.pick(field = "HGT") ... % HGT (height) data field
    .index(time = 2) ... % second time (3-hour forecast, in this case)
    .range(lat = [24 50], lon = [-66 -120]); % lat/lon range of the contiguous U.S.

% create object to query by lat/lon/pressure
height_at_pos = height_data.interpolant(["lat", "lon", "pressure"], "spline");

figure;
surf(height_data.lon, height_data.lat, double(height_data.pick(pressure = 100e3)));
surf(height_data.lon, height_data.lat, double(height_data.pick(pressure = 97.5e3)));
latq = linspace(25, 30, 20);
lonq = linspace(-100, -80, 100);
surf(lonq, latq, double(height_at_pos(latq, lonq, 99e3)));
xlabel("Longitude [deg]");
ylabel("Latitude [deg]");
zlabel("Height [gpm]");
view(3);
profile viewer;

% Exercise: modify above data to color the surface using temperature
% Hint: you need to get the TMP data field
% Note: pay close attention to what data end up as 3-D vs 2-D; if the data has
% multiple fields/pressures/times, you need to select down to a single one for <surf>


%% Math for an air column - Arithmetic and interpolation
% To query along the time axis, it needs to be turned into a numeric axis
% A simple way to do this is to define the time along that axis as 
% "seconds since the start (epoch) of the data"
epoch = min(atmosphere_data.time);
atmosphere_data.time = seconds(atmosphere_data.time - epoch);
atmos_at_point = atmosphere_data.interpolant(["lat", "lon", "time"], "linear");

% Sample the air column at (42.7N 77.1W) interpolated to 1 hour
latq = 42.7;
lonq = -77.1;
forecast = hours(1);
aircolumn = atmos_at_point(latq, lonq, seconds(forecast)); % query interpolant at desired time
aircolumn = squeeze(aircolumn); 
% the output of the interpolant has the dimensions you sampled with, but they
% are scalar values, so get rid of singleton dimensions

% 
H = double(aircolumn.pick(field = "HGT"));
T = aircolumn.pick(field = "TMP");
rho = (1/287.05) * aircolumn.pressure ./ (273.15 + aircolumn.pick(field = "TMP"));
% Notes:
%   - Scalar-vector (e.g. 273 + T) operations always work because the operation
%   doesn't create new dimensions in the xarray
%   - Data axes are always used as column-vectors. The above element-wise
%   vector-vector division only works if aircolumn.pick(field = "TMP") returns a
%   column-vector, which in this case it does, because we used <squeeze>
%   earlier to get rid of the leading scalar dimesions (otherwise it would be
%   1x1x1xN, which is not a compatible size)

figure;
tiledlayout(1,4);

nexttile;
plot(double(T), double(H), "^k");
ylabel("Height [gpm]");
xlabel("Temperature [C]");

nexttile;
plot(aircolumn.pressure, double(H), "^k");
ylabel("Height [gpm]");
xlabel("Pressure [Pa]");

nexttile;
plot(double(rho), double(H), "^k");
ylabel("Height [gpm]");
xlabel("Density [kg/m^3]");

nexttile;
plot(double(aircolumn.pick(field = "UGRD")), double(H), "^", DisplayName = "East wind");
plot(double(aircolumn.pick(field = "VGRD")), double(H), "^", DisplayName = "North wind");

height_sampled = aircolumn.pick(field = ["UGRD", "VGRD"]);
height_sampled = height_sampled.rename("pressure", "height");
height_sampled.height = double(aircolumn.pick(field = "HGT"));
aircolumn_in_height = height_sampled.interpolant("height", "spline");
heightq = linspace(100, 2e3, 200);
fine_sampled = aircolumn_in_height(heightq);

plot(double(fine_sampled.pick(field = "UGRD")), heightq, ...
    "--", DisplayName = "East wind (interpolated)");
plot(double(fine_sampled.pick(field = "VGRD")), heightq, ...
    "--", DisplayName = "North wind (interpolated)");
legend;


