function output_path = try_find_file( search_dir, pattern, msg, allowed_extensions, get_first )
%TRY_FIND_FILE Tries to find a file matching a given pattern in path.
%   Otherwise, prompts the user to select a file matching
%   allowed_extensions, with a dialog box containing msg.

    if nargin < 5
        get_first = false;
        if nargin < 4
            allowed_extensions = '*';
        end
    end

    f=dir(fullfile(search_dir, pattern));

    if get_first && ~isempty(f)
        output_path = [f(1).folder, '\', f(1).name];
    elseif length(f) == 1
        output_path = [f.folder, '\', f.name];
    else
        cwd=pwd();
        try
            cd(search_dir);
        catch
            error([search_dir, ' is not a folder']);
        end

        [fname,outpath] = uigetfile(allowed_extensions, msg);

        if fname == 0
            error('Program terminated. File required to continue.');
        end

        output_path = [outpath,fname];
        cd(cwd);
    end
end

