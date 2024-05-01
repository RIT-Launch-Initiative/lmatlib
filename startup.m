set(groot, "DefaultFigureWindowStyle", "docked");
set(groot, "DefaultAxesFontSize", 14); 
set(groot, "DefaultAxesNextPlot", "add"); % default hold-on

grids = compose("defaultAxes%sGrid", ["X", "Y", "Z"]);
for grid = grids
    set(groot, grid, "on");
end

% replace all interpreters with LaTeX
% taken from MATLAB Answers
%{
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name, 'latex');
end
%}

clear;
