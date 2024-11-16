% Export one figure to a PDF
% figure_to_pdf([name = file_name], [num = figure_number])
%   name:   name of the file to write to (not including .pdf) 
%           Optional - defaults to figure name, errors if figure has no name
%   num:    figure number
%           Optional - defaults to current figure
% Examples:
%   figure_to_pdf
%   figure_to_pdf(num = 2)
%   figure_to_pdf(name = "current figure")
%   figure_to_pdf(num = 1, name = "Figure 1")
function figure_to_pdf(kwargs)
    arguments
        kwargs.name (1,1) string = "";
        kwargs.num (1,1) = -1;
    end
    
    name = kwargs.name;
    num = kwargs.num;

    if num == -1
        fig = gcf();
    else 
        figs = findobj(type = "figure");
        fig = figs([figs.Number] == num);
    end
    if isempty(fig)
        error("Figure %d not found", num);
    end
    if name == ""
        name = fig.Name;
    end

    if isempty(name)
        error("No file name provided or present in figure");
    end

    filename = sprintf("%s.pdf", name);
    fprintf("Saving figure %d to path %s\n", fig.Number, filename);
    exportgraphics(fig, filename, ContentType = "vector", BackgroundColor = "none");
end
