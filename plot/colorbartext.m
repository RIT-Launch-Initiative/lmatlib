function [cb] = colorbartext(label, varargin)
    % Create colorbar with label
    % cb = colorbartext(label, ...)
    cb = colorbar(varargin{:});
    cb.Label.String = label;
end

