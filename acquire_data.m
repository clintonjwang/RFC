function acquire_data(patients, data_dir, working_dir, train_bool, filename_map)
%ACQUIRE_DATA(patients, path, train_bool)
% patients should be cell array of subfolder names each containing a patient
% path is the path to the patient folders
% set train_bool to true if training
        
    disp('Acquiring patient data...');
    parfor i=1:length(patients)
        if ~exist([working_dir,'/data_',patients{i},'.mat'],'file')
            pat = patients{i};
            data_i = acquire_data_single_pat([data_dir,'/',pat], train_bool, filename_map);
            save_wrapper(data_i, [working_dir,'/data_',pat,'.mat'])
        end
    end
end

function data = acquire_data_single_pat(data_dir, train_bool, filename_map)
%acquire_data_single_pat(path, train_bool)
% data_dir is the path to the patient folders
% set train_bool to true if training
% filename_map specifies the file path pattern for each required image file

    R=1.75; %desired low-resolution in mm
    nii_ext = {'*.nii; *.hdr; *.img; *.nii.gz'};
    
    %% Initialize maps
%     niiname_map, isoname_map, longname_map = init_maps(data_dir, train_bool, filename_map);
    niiname_map = containers.Map;
    isoname_map = containers.Map;
    longname_map = containers.Map;
    
    longname_map('pre') = 'T1 pre-contrast image';
    longname_map('art') = 'T1 arterial phase';
    longname_map('pv') = 'T1 portal venous phase';
    longname_map('t2') = 'T2 image';
    longname_map('liver_seg') = 'whole liver segmentation';
    
    for k = {'pre','art','pv','t2'}
        k = char(k);
        if contains(filename_map(k),'.nii')
            niiname_map(k) = try_find_file(data_dir, filename_map(k),...
                ['Select the ', longname_map(k)], nii_ext);
        else
            niiname_map(k) = [data_dir, '/', k, '.nii'];
            fn = dicm2nii(filename_map(k), data_dir, 0);
            movefile([data_dir, '/', fn{1}, '.nii'], niiname_map(k));
%             [V,~,~] = dicomreadVolume(filename_map(k));
%             tmp = double(flip_image(squeeze(V)));
%             save_nii(tmp, niiname_map(k));
        end
    end
    
%     niiname_map('t1_bfc') = try_find_file(data_dir, filename_map('t1-bfc'),...
%                 'Select the T1 bias field estimate', nii_ext);
%     isoname_map('t1_bfc') = '/temp/bias_field_isotropic.nii';
%     isoname_map('t1_bfc') = try_find_file(data_dir, filename_map('t1-bfc'),...
%                 'Select the T1 bias field estimate', nii_ext);
    if isKey(filename_map, 't1_bfc')
        isoname_map('t1_bfc') = filename_map('t1_bfc');
    end
    isoname_map('pre') = '/temp/pre_isotropic.nii';
    isoname_map('art') = '/temp/20s_isotropic.nii';
    isoname_map('pv') = '/temp/70s_isotropic.nii';
    isoname_map('t2') = '/temp/t2_isotropic.nii';
            
    niiname_map('liver_seg') = fullfile(data_dir, '/temp/whole_liver.nii');
    
    isoname_map('liver_seg') = '/temp/whole_liver_isotropic.nii';

    if train_bool
        segs = {'liver_seg', 'vasc_seg', 'tumor_seg', 'necro_seg'};
        
        longname_map('tumor_seg') = 'tumor segmentation';
        longname_map('vasc_seg') = 'vasculature segmentation';
        longname_map('necro_seg') = 'necrosis segmentation';
        
        niiname_map('tumor_seg') = fullfile(data_dir, '/temp/tumor.nii');
        niiname_map('vasc_seg') = fullfile(data_dir, '/temp/vessel.nii');
        niiname_map('necro_seg') = fullfile(data_dir, '/temp/necrosis.nii');
        
        isoname_map('tumor_seg') = '/temp/tumor_isotropic.nii';
        isoname_map('vasc_seg') = '/temp/vessel_isotropic.nii';
        isoname_map('necro_seg') = '/temp/necrosis_isotropic.nii';
    else
        segs = {'liver_seg'};
    end
    
    %% Make isotropic niis
    [~,~,~] = mkdir([data_dir,'/temp']);
    
    % get pixel dims from arterial phase nifti
    art = load_nii(niiname_map('art'));
    temp_res(1) = art.hdr.dime.pixdim(2);
    temp_res(2) = art.hdr.dime.pixdim(3);
    temp_res(3) = art.hdr.dime.pixdim(4);
    [N1,N2,N3] = size(double(flip_image(art.img)));
    data.orig_dims = [N1 N2 N3];

    %make niis from segs
    for seg = 1:length(segs)
        seg = segs{seg};
        if endsWith(data_dir, '/0')
            f=filename_map(seg);
        else
            f=try_find_file(data_dir, filename_map(seg),...
                        ['Select the ', longname_map(seg)], '*.ids');
        end
        mask = get_mask(f, N1,N2,N3);
        nii = make_nii(flip_image(mask),temp_res);
        save_nii(nii,niiname_map(seg));
    end

    % make isotropic niis
    verbose = false;
    for k = keys(niiname_map)
        k = char(k);
        reslice_nii(niiname_map(k),...
            [data_dir,isoname_map(k)],...
            [R,R,R], verbose);
    end
    
    art = load_nii(fullfile(data_dir, isoname_map('art')));
    t2 = load_nii(fullfile(data_dir, isoname_map('t2')));
    if any(size(art.img) ~= size(t2.img))
        [optimizer, metric] = imregconfig('multimodal');
        t2.img = imregister(t2.img, art.img, 'affine', optimizer, metric);
        t2.hdr.dime.dim(2:4) = size(t2.img);
        save_nii(t2, fullfile(data_dir, isoname_map('t2')));
    end
%     end

    %% Load niis into data struct
    for k = {'pre','pv','t2','t1_bfc'}
        k = char(k);
        if isKey(isoname_map,k)
            data.(k) = get_nii_img(fullfile(data_dir, isoname_map(k)));
        end
    end
    for k = 1:length(segs)
        k = segs{k};
        data.(k) = get_nii_img(fullfile(data_dir, isoname_map(k))) > 0;
    end
    data.art = load_nii(fullfile(data_dir, isoname_map('art')));
    data.resolution(1) = data.art.hdr.dime.pixdim(2); %extract voxel dims
    data.resolution(2) = data.art.hdr.dime.pixdim(3);
    data.resolution(3) = data.art.hdr.dime.pixdim(4);
    data.art = double(flip_image(data.art.img));

    %% Crop image to liver mask
    %shrink the masks to proper size if necessary
    [N1,N2,N3] = size(data.art);
    data = shrink_masks(data, N1, N2, N3, train_bool);

    %compute tightened liver mask
    r=4;
    rs=r/data.resolution(1);
    D=bwdist(1-data.liver_seg);
    data.tight_liver_mask = zeros(N1,N2,N3);
    data.tight_liver_mask(D>rs) = data.liver_seg(D>rs);
    
    if train_bool
        data.necro_seg = data.necro_seg .* data.tight_liver_mask;
        data.tumor_seg = data.tumor_seg .* data.tight_liver_mask .*(1-data.necro_seg);
        data.vasc_seg = data.vasc_seg .* data.tight_liver_mask .* (1-data.tumor_seg) .* (1-data.necro_seg);
    end

    %find edges of the liver and crop the images to this dimension
    [i,j,k] = ind2sub(size(data.liver_seg), find(data.liver_seg));
    data.cutoffs_high = [max(i), max(j), max(k)];
    data.cutoffs_low = [min(i), min(j), min(k)];

    for sf={'pre', 'art', 'pv', 't2', 't1_bfc',...
        'tight_liver_mask', 'liver_seg', 'vasc_seg', 'tumor_seg', 'necro_seg'}
        if(isfield(data, char(sf)))
            data.(char(sf)) = data.(char(sf))(data.cutoffs_low(1):data.cutoffs_high(1), ...
                data.cutoffs_low(2):data.cutoffs_high(2),data.cutoffs_low(3):data.cutoffs_high(3));
        end
    end
    
    %% Cleanup
    data.liver_contour = get_contour(data.liver_seg);
    data.tight_liver_contour = get_contour(data.tight_liver_mask);
    if train_bool
        data.tumor_contour = get_contour(data.tumor_seg);
        data.vasc_contour = get_contour(data.vasc_seg);
        data.necro_contour = get_contour(data.necro_seg);
    end

    [~, ~, ~] = rmdir([data_dir,'/temp'], 's');
end

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

end

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

end

function image_o = flip_image(image)

    for k=1:size(image,3)
        image_o(:,:,k) = fliplr(fliplr(image(:,:,k))');
    end

end

function img = get_nii_img(nii_path)
    nii = load_nii(nii_path);
    img = double(flip_image(nii.img));
end

function mask = get_nii_mask(nii_path)
    img = get_nii_img(nii_path);
    mask = img>0;
end

function data = shrink_masks(data, N1, N2, N3, train_bool)
%SHRINK_MASKS Summary of this function goes here
%   Detailed explanation goes here
    if train_bool
        segs = {'liver_seg', 'vasc_seg', 'tumor_seg', 'necro_seg'};
    else
        segs = {'liver_seg'};
    end

    for seg = 1:length(segs)
        seg = segs{seg};
        
        temp = zeros(size(data.art));
        x_max = min(size(data.(seg),1),N1); 
        y_max = min(size(data.(seg),2),N2);
        z_max = min(size(data.(seg),3),N3);
        temp(1:x_max,1:y_max,1:z_max) = data.(seg)(1:x_max,1:y_max,1:z_max);
        data.(seg) = temp;
        if strcmp(seg, 'liver_seg')
            data.iso_full_size = [x_max y_max z_max];
        end
    end
end