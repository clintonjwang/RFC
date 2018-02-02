function restore_mask( nii_path, save_dir, img_name )
%RESTORE_MASK Summary of this function goes here
%   Detailed explanation goes here

    nii = load_nii(nii_path);
    img = double(flip_image(nii.img));
    mask = img>0;
    [mask_root,mask_name,~] = fileparts(nii_path);
    mask_root = fullfile(mask_root, '..');
    
    if nargin == 1
        save_dir = mask_root;
    end
    
    write_ids_mask( mask, mask_root, save_dir, mask_name );
    
    if nargin == 3
%         display_masked_img(img_name, {mask_name}, save_dir);
        display_scrolling_mask(img_name, save_dir, save_dir, {mask_name}, {mask_name})
    end
end