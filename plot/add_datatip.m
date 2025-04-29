%% Add data tip to plot
% add_datatip(plot_handle, label, value[, format])
% Inputs
%   plot_handle: line handle (output argument of plot())
%   label, value, format: see dataTipTextRow documentation
%   
% NOTE: <value> can be scalar; it is internally repeated to a compatible size

function add_datatip(plot_handle, varargin)
    if isscalar(varargin{2})
        varargin{2} = repmat(varargin{2}, size(plot_handle.XData));
    end
    row = dataTipTextRow(varargin{:});
    plot_handle.DataTipTemplate.DataTipRows(end+1) = row;
end
