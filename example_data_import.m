clear; close all;
import lplot.*;

%% RAW DATA IMPORT AND PLOT EXAMPLE
% Pressure to altitude conversion function
H = @(P) (1 - (P / 1013.25) .^ 0.190284) * 145366.45;

% Import data - don't fix variable names with spaces
grim_slow = readtable("sample_data\L1 Grim Test 1Hz.csv", VariableNamingRule = "preserve");
grim_fast = readtable("sample_data\L1 Grim Test 500Hz.csv", VariableNamingRule = "preserve");
rrc3 = readtable("sample_data\L1 RRC3.csv", VariableNamingRule = "preserve");
rs41 = readtable("sample_data\L1 RS41 Test.csv", VariableNamingRule = "preserve");
% MATLAB imports table files as Tables, which you can think of as just structures where each field is a column
% i.e. a column named "pressure" can be accessed as grim_slow.pressure, and will be a column vector

% Create timetables from each table
% field names (.pressure, .temp) are the same as in the CSV

rrc3_t = timetable(rrc3.Pressure, rrc3.Temperature, rrc3.Altitude, ...
    RowTimes = seconds(rrc3.Time));
rrc3_t.Properties.VariableNames = ["Pressure", "Temperature", "Altitude"];
rrc3_t.Properties.VariableUnits = ["mbar", "deg C", "m (AGL)"];

% how much time to take off of Grim's data to line up with RRC3
grim_time_corr = seconds(94.25 - 94.7);
% we only care about pressure and temperature from Grim, so take those directly
grim_slow_t = timetable(grim_slow.pressure, grim_slow.temp, ...
    RowTimes = seconds(grim_slow.timestamp) + grim_time_corr); % seconds(...) tells MATLAB how to interpret the number as a duration
grim_slow_t.Properties.VariableNames = ["Pressure", "Temperature"]; % re-assign the variable names to be nicer
grim_slow_t.Properties.VariableUnits = ["kPa", "deg C"]; % assign units 

% convert kpa to mbar
% the table needs to have units for us to do this
grim_slow_t = convert(grim_slow_t, "kPa", "mbar", @(p) p * 10);
% "@(p) p / 10" is an anonymous function handle for a function that just divides its input by 10
% if you wanted to convert C to F, you could do <convert(rrc3_t, "deg C", "deg F", @(T) T * (5/9) + 32)>

% add column for altitude and assign its units
grim_slow_t.Altitude = H(grim_slow_t.Pressure);
grim_slow_t.Properties.VariableUnits(3) = "ft (ASL)"; % (3) because it's the 3rd column
grim_slow_t = convert(grim_slow_t, "ft (ASL)", "m (AGL)", ...
    @(h) 0.3048 * (h - grim_slow_t.Altitude(end))); % subtract off the ground level (last data point) and convert to m


rs41_time_corr = seconds(90.05 - 88.5);
rs41_t = timetable(rs41.alt, RowTimes = seconds(rs41.Time) + rs41_time_corr);
rs41_t.Properties.VariableNames = "Altitude";
rs41_t.Properties.VariableUnits = "m (ASL)";

rs41_t = convert(rs41_t, "m (ASL)", "m (AGL)", @(h) h - rs41_t.Altitude(1));


%% OpenRocket Import
omen = import_openrocket("sample_data/omen_30mph.csv", verbose = true);

%% Save processed data
% so you don't need to re-import from raw files
save("sample_data/grim_test.mat", "grim_slow_t", "rs41_t", "rrc3_t");
save("sample_data/omen.mat", "omen");
