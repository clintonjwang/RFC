function write_ids_mask( mask, search_path, save_dir, name )
%MAKE_IDS_MASK Summary of this function goes here
%   write_ids_mask( mask, search_path, save_dir, name )

    fpath = try_find_file(search_path, '**/*.ics',...
                'Select any .ics file for this study', '*.ics', true);
    ids_filename = [save_dir, '\', name, '.ids'];
    fileID = fopen(fpath);
    A = fscanf(fileID,'%c');
    ics_text = replace(A, A(regexp(A, '\nfilename')+1:regexp(A, '.ids')+3),...
        [sprintf('filename\t'), pwd(), '\', ids_filename]);

    fileID = fopen([save_dir, '\', name, '.ics'],'w');
    fwrite(fileID, ics_text);

    mask = uint8(transpose_mask_slices(mask, 'w'));
    fileID = fopen(ids_filename,'w');
    fwrite(fileID, mask);
    fclose('all');

end

