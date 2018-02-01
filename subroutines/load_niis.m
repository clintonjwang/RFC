function data = load_niis( data, data_dir, train_bool, use_bias_field, isoname_map )
%LOAD_NIIS Summary of this function goes here
%   Detailed explanation goes here

nii_ext = {'*.nii; *.hdr; *.img; *.nii.gz'};
                
for k in 
data.pre = get_nii_img(isoname_map('pre'));

%load 20s image
data.art = load_nii(try_find_file(data_dir, 'temp/20s_isotropic.nii*',...
            'Select the arterial phase', nii_ext));
data.resolution(1) = data.art.hdr.dime.pixdim(2); %extract voxel dims
data.resolution(2) = data.art.hdr.dime.pixdim(3);
data.resolution(3) = data.art.hdr.dime.pixdim(4);
data.art = double(flip_image(data.art.img));

%load 70s image 
data.pv = get_nii_img(try_find_file(data_dir, 'temp/70s_reg_isotropic.nii',...
            'Select the venous phase', nii_ext));

%load t2 image
data.t2 = get_nii_img(try_find_file(data_dir, 'temp/t2_bfc_reg_isotropic.nii',...
            'Select the T2 image', nii_ext));

%load T1 bias field estimate 
if use_bias_field
    data.bf = get_nii_img(try_find_file(data_dir, 'temp/bias_field_isotropic.nii',...
                'Select the T1 bias field estimate', nii_ext));
end

%load liver_mask
data.liver_mask = get_nii_img(try_find_file(data_dir, 'temp/whole_liver_isotropic.nii',...
            'Select the whole liver mask', nii_ext));
data.liver_mask = data.liver_mask>0;

if train_bool
    %load vessel mask 
    data.vessel_mask = get_nii_img(try_find_file(data_dir, 'temp/vessel_isotropic.nii',...
            'Select the whole liver mask', nii_ext));
    data.vessel_mask = data.vessel_mask>0;

    %load necrosis mask 
    data.necrosis_mask = get_nii_img(try_find_file(data_dir, 'temp/necrosis_isotropic.nii',...
            'Select the necrosis mask', nii_ext));
    data.necrosis_mask = data.necrosis_mask>0;

    %load tumor mask 
    data.tumor_mask = get_nii_img(try_find_file(data_dir, 'temp/tumor_isotropic.nii',...
            'Select the tumor mask', nii_ext));
    data.tumor_mask = data.tumor_mask>0;
end

return

end

function img = get_nii_img(nii_path)
    nii = load_nii(nii_path);
    img = double(flip_image(nii.img));
end
