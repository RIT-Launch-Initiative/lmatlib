%% Automatically export all figures to PDFs
% figures_to_pdfs()
%   Uses the name of the figure as the file name
%   All figures have to have names and they must be unique
function figures_to_pdfs
    figs = findobj(type = "figure")';
    names = {figs.Name};

    if any(names == "")
        error("All figures must have names");
    end
    if numel(names) ~= numel(unique(names))
        error("Figure names must be unique");
    end

    for fig = figs
        filename = sprintf("%s.pdf", fig.Name);
        fprintf("Saving figure %d to path %s\n", fig.Number, filename);
        exportgraphics(fig, filename, ContentType = "vector", BackgroundColor = "none");
    end
end
