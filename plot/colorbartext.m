function [cb] = colorbartext(label, varargin)
    cb = colorbar(varargin{:});
    cb.Label.String = label;
end

