function data = acquire_data_single_pat(data_dir, train_bool, use_bias_field)
%acquire_data_single_pat(path, train_bool)
% path is the path to the patient folders
% set train_bool to true if training

R=1.75; %desired low-resolution in mm

nii_ext = {'*.nii; *.hdr; *.img; *.nii.gz'};

data.pre = load_nii(try_find_file(data_dir, '**/pre.nii',...
                    'Select the pre-contrast nifti file', nii_ext));

mkdir([data_dir,'/temp']);
% fpath=dir(fullfile([path, patient, segs_dir], '**/wholeliver.ids'));
if true%exist([data_dir,'/temp/whole_liver.nii'],'file') == 0
    % get pixel dims from arterial phase nifti
    temp = load_nii(try_find_file(data_dir, '**/20s.nii',...
                    'Select the arterial phase nifti file', nii_ext));
    temp_res(1) = temp.hdr.dime.pixdim(2); %extract pixel volume
    temp_res(2) = temp.hdr.dime.pixdim(3);
    temp_res(3) = temp.hdr.dime.pixdim(4);

    [N1,N2,N3] = size(double(flip_image(temp.img)));
    data.orig_dims = [N1 N2 N3];

    %make nii file from whole liver segmentation
    f=try_find_file(data_dir, '**/whole*liver.ids',...
                'Select the whole liver segmentation', '*.ids');
    liver_mask = get_mask(f, N1,N2,N3);
    liver_nii = make_nii(flip_image(liver_mask),temp_res);
    save_nii(liver_nii,[data_dir,'/temp/whole_liver.nii']);

    if train_bool
        %flags indicating the existence of vessel, tumor, and necrosis
        %segmentations 
        f=try_find_file(data_dir, '**/*vessel*.ids',...
                    'Select the vasculature segmentation', '*.ids');
        mask = get_mask(f, N1,N2,N3);
        nii = make_nii(flip_image(mask),temp_res);
        save_nii(nii,[data_dir,'/temp/vessel.nii']);
        
        f=try_find_file(data_dir, '**/*tumor*.ids',...
                    'Select the tumor segmentation', '*.ids');
        mask = get_mask(f, N1,N2,N3);
        nii = make_nii(flip_image(mask),temp_res);
        save_nii(nii,[data_dir,'/temp/tumor.nii']);
        
        f=try_find_file(data_dir, '**/*nec*.ids',...
                    'Select the tumor segmentation', '*.ids');
        mask = get_mask(f, N1,N2,N3);
        nii = make_nii(flip_image(mask),temp_res);
        save_nii(nii,[data_dir,'/temp/necrosis.nii']);
    end

    make_isotropic_niis(data_dir, R, train_bool, use_bias_field);
end

data = load_niis(data, data_dir, train_bool, use_bias_field);

%get T1 image dimensions
[N1,N2,N3] = size(data.art);

%shrink the masks to proper size if necessary
data = shrink_masks(data, N1, N2, N3, train_bool);

%compute tightened liver mask
r=4;
rs=r/data.resolution(1);
D=bwdist(1-data.liver_mask);
data.tight_liver_mask = zeros(N1,N2,N3);
data.tight_liver_mask(D>rs) = data.liver_mask(D>rs);

if train_bool
    data.necrosis_mask = data.necrosis_mask .* data.tight_liver_mask;
    data.tumor_mask = data.tumor_mask .* data.tight_liver_mask .*(1-data.necrosis_mask);
    data.vessel_mask = data.vessel_mask.*data.tight_liver_mask.*(1-data.tumor_mask).*(1-data.necrosis_mask);
end
    
%find edges of the liver and crop the images to this dimension
[i,j,k] = ind2sub(size(data.liver_mask),find(data.liver_mask));
i_min = min(i);
i_max = max(i);

j_min=min(j);
j_max=max(j);

k_min=min(k);
k_max=max(k);

data.cutoffs_high = [i_max, j_max, k_max];
data.cutoffs_low = [i_min, j_min, k_min];

sf={'pre','art','pv','t2','liver_mask','tumor_mask','bf',...
    'tight_liver_mask','vessel_mask','necrosis_mask'};

for sf_count=1:length(sf)
    if(isfield(data,sf{sf_count}))
        data.(sf{sf_count}) = data.(sf{sf_count})(i_min:i_max, ...
            j_min:j_max,k_min:k_max);
    end
end

%get contours for later visualization
data.liver_contour = get_contour(data.liver_mask);
data.tight_liver_contour = get_contour(data.tight_liver_mask);
if train_bool
    data.tumor_contour = get_contour(data.tumor_mask);
    data.vessel_contour = get_contour(data.vessel_mask);
    data.necrosis_contour = get_contour(data.necrosis_mask);
end

[~, ~, ~] = rmdir([data_dir,'/temp'], 's');

return

function [mask] = get_mask(fpath,N1,N2,N3)

fileID = fopen(fpath);
A = fread(fileID);
ind = find(A);
[i, j, k]=ind2sub([N2,N1,N3],ind);
fclose('all');

mask = zeros(N1,N2,N3);
for count=1:length(i)
    mask(i(count),j(count),k(count))=1;
end

mask = transpose_mask_slices(mask, 'r');

return

function [contour] = get_contour(mask)

[N1,N2,N3] = size(mask);

for n=1:N3
    perim(:,:,n) = bwperim(mask(:,:,n));
end

for n3=1:N3
    count=1;
    coordinates=[];
    for n1=1:N1
        for n2=1:N2
            if(perim(n1,n2,n3)==1)
                coordinates(count,1)=n1;
                coordinates(count,2)=n2;
                count = count + 1;
            end
        end
    end
    contour{n3} = coordinates;
end

return

function data = shrink_masks(data, N1, N2, N3, train_bool)
    
temp = zeros(size(data.art)); 
x_max=min(size(data.liver_mask,1),N1); 
y_max=min(size(data.liver_mask,2),N2);
z_max=min(size(data.liver_mask,3),N3);
temp(1:x_max,1:y_max,1:z_max)=data.liver_mask(1:x_max,1:y_max,1:z_max);
data.iso_full_size = [x_max y_max z_max];
data.liver_mask=temp;

if train_bool
    temp = zeros(size(data.art)); 
    x_max=min(size(data.vessel_mask,1),N1); 
    y_max=min(size(data.vessel_mask,2),N2);
    z_max=min(size(data.vessel_mask,3),N3);
    temp(1:x_max,1:y_max,1:z_max)=data.vessel_mask(1:x_max,1:y_max,1:z_max);  
    data.vessel_mask=temp;

    temp = zeros(size(data.art)); 
    x_max=min(size(data.necrosis_mask,1),N1); 
    y_max=min(size(data.necrosis_mask,2),N2);
    z_max=min(size(data.necrosis_mask,3),N3);
    temp(1:x_max,1:y_max,1:z_max)=data.necrosis_mask(1:x_max,1:y_max,1:z_max);  
    data.necrosis_mask=temp;

    temp = zeros(size(data.art)); 
    x_max=min(size(data.tumor_mask,1),N1); 
    y_max=min(size(data.tumor_mask,2),N2);
    z_max=min(size(data.tumor_mask,3),N3);
    temp(1:x_max,1:y_max,1:z_max)=data.tumor_mask(1:x_max,1:y_max,1:z_max);  
    data.tumor_mask=temp;
end

return

function image_o = flip_image(image)%

for k=1:size(image,3)
    image_o(:,:,k) = fliplr(fliplr(image(:,:,k))');
end

return