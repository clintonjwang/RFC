function data = load_niis( data, data_dir, train_bool, use_bias_field )
%LOAD_NIIS Summary of this function goes here
%   Detailed explanation goes here

nii_ext = {'*.nii; *.hdr; *.img; *.nii.gz'};
                
%load pre image 
data.pre = load_nii(try_find_file(data_dir, '**/pre_reg_isotropic.nii',...
            'Select the pre-contrast image', nii_ext));
data.pre = double(flip_image(data.pre.img));

%load 20s image
data.art = load_nii(try_find_file(data_dir, '**/20s_isotropic.nii*',...
            'Select the arterial phase', nii_ext));
data.resolution(1) = data.art.hdr.dime.pixdim(2); %extract voxel dims
data.resolution(2) = data.art.hdr.dime.pixdim(3);
data.resolution(3) = data.art.hdr.dime.pixdim(4);
data.art = double(flip_image(data.art.img));

%load 70s image 
data.pv = load_nii(try_find_file(data_dir, '**/70s_reg_isotropic.nii',...
            'Select the venous phase', nii_ext));
data.pv = double(flip_image(data.pv.img));

%load t2 image
data.t2 = load_nii(try_find_file(data_dir, '**/t2_bfc_reg_isotropic.nii',...
            'Select the T2 image', nii_ext));
data.t2 = double(flip_image(data.t2.img));

%load T1 bias field estimate 
if use_bias_field
    data.bf = load_nii(try_find_file(data_dir, '**/bias_field_isotropic.nii',...
                'Select the T1 bias field estimate', nii_ext));
    data.bf = double(flip_image(data.bf.img));
end

%load liver_mask
data.liver_mask = load_nii(try_find_file(data_dir, '**/whole_liver_isotropic.nii',...
            'Select the whole liver mask', nii_ext));
data.liver_mask = double(flip_image(data.liver_mask.img));
data.liver_mask = data.liver_mask>0;

if train_bool
    %load vessel mask 
    data.vessel_mask = load_nii(try_find_file(data_dir, '**/vessel_isotropic.nii',...
            'Select the whole liver mask', nii_ext));
    data.vessel_mask = double(flip_image(data.vessel_mask.img));
    data.vessel_mask = data.vessel_mask>0;

    %load necrosis mask 
    data.necrosis_mask = load_nii(try_find_file(data_dir, '**/necrosis_isotropic.nii',...
            'Select the necrosis mask', nii_ext));
    data.necrosis_mask = double(flip_image(data.necrosis_mask.img));
    data.necrosis_mask = data.necrosis_mask>0;

    %load tumor mask 
    data.tumor_mask = load_nii(try_find_file(data_dir, '**/tumor_isotropic.nii',...
            'Select the tumor mask', nii_ext));
    data.tumor_mask = double(flip_image(data.tumor_mask.img));
    data.tumor_mask = data.tumor_mask>0;
end

return

end

