function multiplicity = check_for_file( search_dir, pattern )
%CHECK_FOR_FILE Tries to find a file matching a given pattern in path.
%   If it exists, returns true. Otherwise, false.

    f=dir(fullfile(search_dir, pattern));

    if length(f) == 1
        multiplicity = 1;
    elseif ~isempty(f)
        multiplicity = 2; %[f(1).folder, '\', f(1).name]
    else
        multiplicity = 0;
    end
end
