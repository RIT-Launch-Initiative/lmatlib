%% Get the path of a file relative to the current project
% out = pfullfile(folder, subfolder, ....)
% Read as (project)fullfile
% Wraps fullfile() to pre-prend the current project's root folder to file name generation
function [out] = pfullfile(varargin)
    proj = matlab.project.currentProject;
    out = fullfile(proj.RootFolder, varargin{:});
end

