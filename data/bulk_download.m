% Bulk download files, checking for existing ones, with informative display output
% [files, notfiles] = bulk_download(dest, urls, names[, reader])
% INPUTS
%   dest        destination folder
%   urls        source URLs
%   names       destination file names - appended to <dest>
%   reader      function handle (name, url) -> logical (defaults to websave)
% OUTPUTS
%   files       files successfully downloaded
%   notfiles    files not successfully downloaded

function [files, notfiles] = bulk_download(dest, urls, names, reader)
    arguments
        dest (1, 1) string {mustBeFolder};
        urls (:, 1) string;
        names (:, 1) string;
        reader (1, 1) function_handle = @websave;
    end 

    if length(urls) ~= length(names)
        error("Each URL must have a name");
    end

    files = fullfile(dest, names);
    files = strrep(files, "\", "/");

    % check presence of every file
    on_disk = false(size(files));
    on_server = false(size(files));
    sizes = zeros(size(files));
    for i_file = 1:length(files)
        info = dir(files(i_file));
        if isscalar(info)
            on_disk(i_file) = true;
            sizes(i_file) = info.bytes;
            continue;
        end
        
        info = dir(urls(i_file));
        if isscalar(info)
            on_server(i_file) = true;
            sizes(i_file) = info.bytes;
            continue;
        end
    end

    % emit status messages for files located on disk, on server, or not at all
    n_on_disk = sum(on_disk);
    if n_on_disk > 0
        fprintf("%d files located on disk:\n", n_on_disk);
        fprintf(listing(files(on_disk), sizes(on_disk)));
    end

    n_on_server = sum(on_server);
    if n_on_server > 0
        fprintf("%d files located on server:\n", n_on_server);
        fprintf(listing(urls(on_server), sizes(on_server)));
    end

    n_nowhere = sum(~on_disk & ~on_server);
    if n_nowhere > 0
        files_nowhere = urls(~on_disk & ~on_server);
        fprintf(2, "%d files not found:\n%s", n_nowhere, compose("\t%s\n", files_nowhere).join(""));
    end

    % download the files located on server
    i_on_server = find(on_server); % linear indices of files on server
    for i_file = 1:length(i_on_server) % download the file described by each linear index
        n_file = i_on_server(i_file);
        file_url = urls(n_file);
        fprintf("%d of %d: %s from %s... ", i_file, length(i_on_server), ...
            format_filesize(sizes(n_file)), file_url)

        try
            tic;
            downloaded = reader(files(n_file), file_url);
            timed = toc;
            if downloaded % print status if done or not
                on_disk(n_file) = true;
                fprintf("succeeded in %.1f sec\n", timed);
            else
                fprintf("failed\n");
            end
        catch mex % don't stop for anything
            fprintf("failed\n");
            warning(getReport(mex));
        end

    end

    notfiles = files(~on_disk);
    files = files(on_disk);

    % Helper functions

    % engineering format for file sizes
    function str = format_filesize(sizes)
        pwr = floor(log10(sizes)/3)*3;
        units = dictionary([-Inf 0 3 6 9 12], ["B ", "B ", "KB", "MB", "GB", "TB"]);
        reduced_sizes = sizes./10.^pwr;
        reduced_sizes(~isfinite(reduced_sizes)) = 0;
        str = compose("%5.1f %s", reduced_sizes, units(pwr));
    end

    % create directory listing
    function [str] = listing(files, sizes)
        assert(length(files) == length(sizes));

        [paths, names, exts] = fileparts(files);
        names = names + exts;
        folders = unique(paths);

        fmt = "\t(%s)\t%s";
        if isscalar(folders)
            str = [];
        else
            str = sprintf(fmt, format_filesize(sum(sizes)), ...
                sprintf("TOTAL %d FILES", length(files)));
        end
        for i_folder = 1:length(folders)
            folder = folders(i_folder);
            names_under = names(paths == folder);
            sizes_under = sizes(paths == folder);
            str = join([str;
                sprintf(fmt, format_filesize(sum(sizes_under)), ...
                sprintf("%d files under %s", sum(paths == folder), folder + "/"));
                compose(fmt, format_filesize(sizes_under), names_under)], "\n");
        end
        str = str + "\n";
    end
end
