function data = acquire_data_single_pat(data_dir, train_bool, filename_map)
%acquire_data_single_pat(path, train_bool)
% path is the path to the patient folders
% set train_bool to true if training

    R=1.75; %desired low-resolution in mm

    nii_ext = {'*.nii; *.hdr; *.img; *.nii.gz'};
    niiname_map, isoname_map, longname_map = init_maps(data_dir, train_bool, filename_map);
    niiname_map = containers.Map;
    isoname_map = containers.Map;
    longname_map = containers.Map;
    
    niiname_map('pre') = try_find_file(data_dir, filename_map('pre'),...
                        'Select the pre-contrast nifti file', nii_ext);
    niiname_map('art') = try_find_file(data_dir, filename_map('art'),...
                        'Select the arterial phase nifti file', nii_ext);
    niiname_map('liver_seg') = '/temp/whole_liver.nii';
    
    isoname_map('pre') = '/temp/pre_reg_isotropic.nii';
    isoname_map('art') = '/temp/20s_isotropic.nii';
    isoname_map('ven') = '/temp/70s_reg_isotropic.nii';
    isoname_map('t2') = '/temp/t2_bfc_reg_isotropic.nii';
    isoname_map('liver_seg') = '/temp/whole_liver_isotropic.nii';
    isoname_map('t1_bfc') = try_find_file(data_dir, '**/bias_field_isotropic.nii',...
                'Select the T1 bias field estimate', nii_ext);
    
    longname_map('pre') = 'T1 pre-contrast image';
    longname_map('art') = 'T1 arterial phase';
    longname_map('ven') = 'T1 portal venous phase';
    longname_map('t2') = 'T2 image';
    longname_map('liver_mask') = 'whole liver segmentation';

    if train_bool
        segs = {'liver_seg', 'vasc_seg', 'tumor_seg', 'necro_seg'};
        
        niiname_map('tumor_seg') = '/temp/tumor.nii';
        niiname_map('vasc_seg') = '/temp/vessel.nii';
        niiname_map('necro_seg') = '/temp/necrosis.nii';
        
        isoname_map('tumor_seg') = '/temp/tumor_isotropic.nii';
        isoname_map('vasc_seg') = '/temp/vessel_isotropic.nii';
        isoname_map('necro_seg') = '/temp/necrosis_isotropic.nii';
        
        longname_map('tumor_seg') = 'tumor segmentation';
        longname_map('vasc_seg') = 'vasculature segmentation';
        longname_map('necro_seg') = 'necrosis segmentation';
    else
        segs = {'liver_seg'};
    end
    
    
    mkdir([data_dir,'/temp']);
    
    % get pixel dims from arterial phase nifti
    temp = load_nii(niiname_map('art'));
    temp_res(1) = temp.hdr.dime.pixdim(2);
    temp_res(2) = temp.hdr.dime.pixdim(3);
    temp_res(3) = temp.hdr.dime.pixdim(4);
    [N1,N2,N3] = size(double(flip_image(temp.img)));
    data.orig_dims = [N1 N2 N3];
    
    %make niis from other segs
    for seg = segs
        f=try_find_file(data_dir, filename_map(seg),...
                    ['Select the ', longname_map(seg)], '*.ids');
        mask = get_mask(f, N1,N2,N3);
        nii = make_nii(flip_image(mask),temp_res);
        save_nii(nii,[data_dir,niiname_map(seg)]);
    end

    % make isotropic niis
    verbose = false;
    for k = keys(niiname_map)
        reslice_nii(niiname_map(k),...
            [data_dir,isoname_map(k)],...
            [R,R,R], verbose);
    end
%     end

%     data = load_niis(data, data_dir, train_bool, isoname_map);
    for k = {'pre','ven','t2','t1_bfc'}
        data = setfield(data, get_nii_img(isoname_map(k)));
    end
    for k = segs
        data = setfield(data, get_nii_img(isoname_map(k)));
    end

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

function image_o = flip_image(image)%

for k=1:size(image,3)
    image_o(:,:,k) = fliplr(fliplr(image(:,:,k))');
end

return