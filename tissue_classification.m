function masks = tissue_classification(patients, model_dir, working_dir, train_bool, data_dir, mask_dir, config)
%TISSUE_CLASSIFICATION
%   Detailed explanation goes here
    
    %% Parameters
    disp('Initializing random forest...');
    if train_bool
        save([model_dir, '/tree_config.mat'], 'config');
    else
        load([model_dir, '/tree_config.mat']);
    end

    ntrees = config.ntrees; %trees in each random forest, 800 by default
    RAC = config.RAC; %number of iterations of classification, 3 by default
    sl = config.sl; %structured context patch size, 8 by default
    sl_spherical = config.sl_spherical; %spherical context patch size, 5 by default
    num_bins = config.num_bins; %number of histogram bins for spherical histogram context features, 5 by default
    sa=config.sa; %number of auto-context features to sample per dimension, 6 by default
    min_leaf_size = config.min_leaf_size; %random forest parameter: mininum leaf size, 50 by default
    T_par = 100000; %200000 %number of training samples from parenchyma voxels
    T_other = 30000; %50000 %number of training samples from each non-parenchyma class
    
    num_patients = length(patients); %total number of patients
    num_classes = 4; %number of tissue classes

    %initialize global variables
%     opt = statset('UseParallel',true); %random forest parameter: train trees in parallel 
%     parpool(4); %number of cores to distribute processing over
    M = generate_spherical_masks(sl_spherical); %generate spherical masks
    num_context_features = (sa^3) * (num_classes - 1)... %structured context features
             + ((num_classes-1)*num_bins*length(M)); %plus spherical context features 
    tmp = load_wrapper([working_dir,'/small_features_',patients{1},'.mat']);
    num_intensity_features = tmp.num_intensity_features;
    clear tmp

    voxel_data = load_wrapper([working_dir,'/voxel_data.mat']);
    
    start_r = 3;
    for td=1:num_patients
        if ~exist([working_dir, '/context_',patients{td},'_1.mat'], 'file')
            start_r = 1;
            break
        end
    end
    if start_r > 2
        for td=1:num_patients
            if ~exist([working_dir, '/context_',patients{td},'_2.mat'], 'file')
                start_r = 2;
                break
            end
        end
    end
    
    %% Loop through rounds of auto-context
    for r=start_r:RAC
        %% Train random forest
        if train_bool && ~exist([model_dir,'/tree_',num2str(r),'.mat'], 'file')
            disp('Creating training data samples for random forest...');
            if(r==1)
                num_predictors = num_intensity_features;
            else
                num_predictors = num_intensity_features + num_context_features;
            end
            voxel_data_subset = get_data_subset(T_par, T_other, patients,...
                        voxel_data, r, num_intensity_features, num_context_features, working_dir);
            
            disp(['Fitting random forest (round ', num2str(r), ' of ', num2str(RAC), ')...']);
            for t=1:ntrees
                %bag training data
                num_data_points = round(length(voxel_data_subset.labels) / 3);
                bag = randi(length(voxel_data_subset.labels), [1,num_data_points]);
                if(r==1)
                    tree_data = voxel_data_subset.features(bag,:);
                else
                    tree_data = [voxel_data_subset.features(bag,:), voxel_data_subset.acf(bag,:)];
                end

                C{t} = compact(fitctree(tree_data,voxel_data_subset.labels(bag,:),...
                    'MinLeafSize',min_leaf_size,...
                    'NumVariablesToSample',sqrt(num_predictors),...
                    'SplitCriterion','deviance','Surrogate','on')); %opt
            end
            clear voxel_data_subset;
            save_wrapper(C, [model_dir,'/tree_',num2str(r),'.mat']);
            toc
        else
            disp(['Loading pre-trained forest for round ', num2str(r), '...']);
            C = load_wrapper([model_dir,'/tree_',num2str(r),'.mat']);
        end
        
        %% Classification - compute label probability maps for each tissue class and patient
        disp('Computing tissue classification...');
        locations = load_wrapper([working_dir,'/locations.mat']);
        for td = 1:num_patients %parfor - requires ~15GB RAM to parallelize
            disp(['	', num2str(td), ' of ', num2str(num_patients), '...']);
            f = load_wrapper([working_dir,'/features_',patients{td},'.mat']);
            cls = zeros(numel(f.labels),num_classes);
            if r==1 %if first round of classification: appearance features only 
                for t=1:ntrees
                    [~,scores] = predict(C{t}, f.intensities);
                    cls = cls + scores;
                end
            else %for subsequent rounds, incorporate auto-context features 
                acf = load_context(patients{td}, r-1, num_context_features, working_dir);
                for t=1:ntrees
                    [~,scores] = predict(C{t},[f.intensities, acf]);
                    cls = cls + scores;
                end
            end
%             classification{td} = cls./ntrees;
            cls = cls./ntrees;
            
            if r<RAC
                data = load_wrapper([working_dir,'/norm_data_',patients{td},'.mat']);
                maps = generate_maps_dummy(cls,f.sz,locations{td},num_classes);
                acf = [compute_patches(maps,locations{td},sl,data,sa),...
                    compute_patches_spherical(maps,locations{td},sl_spherical,data,num_bins)];
                save_wrapper(acf, [working_dir, '/context_',patients{td},'_',num2str(r),'.mat']);
            else
                classification{td} = cls;
            end
        end
        toc
        clear C f acf cls locations

        % Compute auto-context features, unless it is the final iteration
%         if r<RAC
%             compute_autocontext_features(classification, sl, sa, patients,...
%                 sl_spherical,num_bins,working_dir,r);
%         end
    end
    
    %% Produce masks
    if ~train_bool
        disp('Creating masks...');
        parfor td=1:num_patients
            f = load_wrapper([working_dir,'/small_features_',patients{td},'.mat']);
            data = load_wrapper([working_dir,'/norm_data_',patients{td},'.mat']);

            scores = classification{td};
            pred_labels = (scores(:,2)>scores(:,1)) .* (scores(:,2)>scores(:,3)).* (scores(:,2)>scores(:,4)) + ...
                 2 * (scores(:,3)>scores(:,1)) .* (scores(:,3)>scores(:,2)).* (scores(:,3)>scores(:,4)) + ...
                 3 * (scores(:,4)>scores(:,1)) .* (scores(:,4)>scores(:,2)).* (scores(:,4)>scores(:,3));
            pred_labels = pred_labels + 1;

            pred_img = zeros(f.sz);
            for pix_idx = 1:length(f.locations)
                pred_img(f.locations(pix_idx)) = pred_labels(pix_idx);
            end

            pad_img = padarray(pred_img, data.iso_full_size - data.cutoffs_high, 0, 'post');
            pad_img = padarray(pad_img, data.cutoffs_low - 1, 0, 'pre');
            rescaled_img = round(imresize3(pad_img, data.orig_dims));
%                 rescaled_img = pred_img;
            nec_mask = rescaled_img == 2;
            vasc_mask = rescaled_img == 3;
            viatumor_mask = rescaled_img == 4;
            paren_mask = rescaled_img == 1;

            masks{td} = {viatumor_mask, nec_mask, vasc_mask, paren_mask};

            [~,~,~]=mkdir(mask_dir);
    %             [vasc_mask, nec_mask, viatumor_mask, paren_mask] = get_masks(classification{td});
            write_ids_mask(vasc_mask, data_dir, mask_dir, 'vasculature_mask');
            write_ids_mask(nec_mask, data_dir, mask_dir, 'necrosis_mask');
            write_ids_mask(viatumor_mask, data_dir, mask_dir, 'viable_tumor_mask');
            write_ids_mask(paren_mask, data_dir, mask_dir, 'parenchyma_mask');
        end
    end
end

function voxel_data_subset = get_data_subset(T_par, T_other, patients, voxel_data, r, nf, nacf, working_dir)
    %get subset of voxel features for this round
    if T_par < length(voxel_data{1}.locs)
        p = randsample(length(voxel_data{1}.locs),T_par);
        normal_locs = voxel_data{1}.locs(p);
        normal_patients = voxel_data{1}.patients(p);
    else
        normal_locs = voxel_data{1}.locs;
        normal_patients = voxel_data{1}.patients;
    end

    if T_other < length(voxel_data{2}.locs)
        p = randsample(length(voxel_data{2}.locs),T_other);
        cancer_locs = voxel_data{2}.locs(p);
        cancer_patients = voxel_data{2}.patients(p);
    else
        cancer_locs = voxel_data{2}.locs;
        cancer_patients = voxel_data{2}.patients;
    end

    if T_other < length(voxel_data{3}.locs)
        p = randsample(length(voxel_data{3}.locs),T_other);
        vessel_locs = voxel_data{3}.locs(p);
        vessel_patients = voxel_data{3}.patients(p);
    else
        vessel_locs = voxel_data{3}.locs;
        vessel_patients = voxel_data{3}.patients;
    end

    if T_other < length(voxel_data{4}.locs)
        p = randsample(length(voxel_data{4}.locs),T_other);
        necrosis_locs = voxel_data{4}.locs(p);
        necrosis_patients = voxel_data{4}.patients(p);
    else
        necrosis_locs = voxel_data{4}.locs;
        necrosis_patients = voxel_data{4}.patients;
    end

    voxel_data_subset.labels = uint8([0*ones(length(normal_locs),1);...
        1*ones(length(cancer_locs),1);2*ones(length(vessel_locs),1);3*ones(length(necrosis_locs),1)]);
    voxel_data_subset.patients = [normal_patients;cancer_patients;vessel_patients;necrosis_patients];
    voxel_data_subset.locs = [normal_locs;cancer_locs;vessel_locs;necrosis_locs];
    voxel_data_subset.features = zeros(length(voxel_data_subset.labels),nf);
    if r > 1
        voxel_data_subset.acf = zeros(length(voxel_data_subset.labels),nacf);
    end
    
    for i=1:length(patients)
        patient_locs = find(voxel_data_subset.patients==i);

        if(~isempty(patient_locs))
            f = load_wrapper([working_dir,'/features_',patients{i},'.mat']);
            for li=patient_locs'
                ind = voxel_data_subset.locs(li);
                voxel_data_subset.features(li,:) = f.intensities(ind,:);
            end
        
            if r > 1
                acf = load_context(patients{i}, r-1, nacf, working_dir);
                for li=patient_locs'
                    ind = voxel_data_subset.locs(li);
                    voxel_data_subset.acf(li,:) = acf(ind,:);
                end
            end
        end
    end
    voxel_data_subset.features = single(voxel_data_subset.features);
    if r > 1
        voxel_data_subset.acf = single(voxel_data_subset.acf);
    end
    toc
end

%dummy function to allow binary file saves inside parfor-loop 
function [] = save_dummy_binary(data,path)
    fileID = fopen(path,'w');
    fwrite(fileID,data,'double');
    fclose(fileID);
end

function [AUC,sensitivity,specificity,precision,accuracy,DSC]...
    = compute_effectivness(scores,labels)
% compute accuracy measures of the classification 

    %loop through labels
    for c=1:3
        %check if patient has corresponding tissue class
        if(~isempty(find(labels==c, 1)))
            lab = double(labels==c);
            [~,~,~,AUC(c,1)] = perfcurve(lab,scores(:,c+1),1);

            if(c==1)
                classification = (scores(:,2)>scores(:,1)) .* (scores(:,2)>scores(:,3)).* (scores(:,2)>scores(:,4));
            elseif(c==2)
                classification = (scores(:,3)>scores(:,1)) .* (scores(:,3)>scores(:,2)).* (scores(:,3)>scores(:,4));
            else
                classification = (scores(:,4)>scores(:,1)) .* (scores(:,4)>scores(:,2)).* (scores(:,4)>scores(:,3));
            end

            precision(c,1)=compute_precision(classification,lab);
            accuracy(c,1)=compute_accuracy(classification,lab);
            sensitivity(c,1)=compute_sensitivity(classification,lab);
            specificity(c,1)=compute_specificity(classification,lab);
            DSC(c,1)=compute_DSC(classification,lab);

        %if not, assign dummy value of -1
        else
            AUC(c,1)=-1;
            precision(c,1)=-1;
            accuracy(c,1)=-1;
            sensitivity(c,1)=-1;
            specificity(c,1)=-1;
            DSC(c,1)=-1;
        end
    end
end

function compute_acf(cls, td, sl, sa, patients, sl_spherical, num_bins, working_dir, r)
    num_classes = 4;
    locations = load_wrapper([working_dir,'/locations.mat']);
    f = load_wrapper([working_dir,'/small_features_',patients{td},'.mat']);
    data = load_wrapper([working_dir,'/norm_data_',patients{td},'.mat']);

    %compute label context features
    maps = generate_maps_dummy(cls,f.sz,locations{td},num_classes);

    acf = [compute_patches(maps,locations{td},sl,data,sa),...
        compute_patches_spherical(maps,locations{td},sl_spherical,data,num_bins)];

    %save auto-context features in binary file 
    save_wrapper(acf, [working_dir, '/context_',patients{td},'_',num2str(r),'.mat']);
%         save_dummy_binary(acf, [working_dir, '/context_',patients{td},'_',num2str(r),'.bin']);
end

function compute_autocontext_features(classification,...
    sl, sa, patients, sl_spherical, num_bins, working_dir, r)

    disp('Extracting auto-context features...')
    num_classes = 4;
    locations = load_wrapper([working_dir,'/locations.mat']);
    
    parfor td=1:length(patients)
        f = load_wrapper([working_dir,'/small_features_',patients{td},'.mat']);
        data = load_wrapper([working_dir,'/norm_data_',patients{td},'.mat']);
        
        %compute label context features
        maps = generate_maps_dummy(classification{td},f.sz,locations{td},num_classes);

        acf = [compute_patches(maps,locations{td},sl,data,sa),...
            compute_patches_spherical(maps,locations{td},sl_spherical,data,num_bins)];

        %save auto-context features in binary file 
        save_wrapper(acf, [working_dir, '/context_',patients{td},'_',num2str(r),'.mat']);
%         save_dummy_binary(acf, [working_dir, '/context_',patients{td},'_',num2str(r),'.bin']);
    end
    toc
end

function acf = load_context(patient, r, num_context_features, working_dir)
    acf = load_wrapper([working_dir, '/context_',patient,'_',num2str(r),'.mat']);
%     f = load_wrapper([working_dir,'/small_features_',patient,'.mat']);
%     fileID = fopen([working_dir, '/context_',patient,'_',num2str(r),'.bin'],'r');
%     acf = fread(fileID,[numel(f.labels), num_context_features],'double');
%     fclose(fileID);
end

%compute structured auto-context features 
function temp = compute_patches(maps,locations,sl,data_i,sa)
    %initialize variables...
    sz=size(data_i.tight_liver_mask);
    [ii,jj,kk] = ind2sub(sz,locations);

    for c=1:size(maps,1)
        for c1=1:size(maps,2)
            patches{c,c1}=zeros(length(ii),sa(c1)^3);
        end
    end
    %loop through data points
    for d=1:length(ii)
        for c1=1:size(maps,2)

            x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,sz(1)),sa(c1)));
            y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,sz(2)),sa(c1)));
            z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,sz(3)),sa(c1)));
            sample_locs = combvec(x,y,z);

            inds = sub2ind(sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));

            for c=1:size(maps,1)
                patches{c,c1}(d,1:(sa(c1))^3) = maps{c,c1}(inds);
            end

        end
    end

    temp=[];
    for c=1:size(maps,1)
        for c1=1:size(maps,2)
            temp=[temp,patches{c,c1}];
        end
    end

    return
end

%compute spherical histogram auto-context features 
function temp = compute_patches_spherical(maps,locations,sl,data_i,num_bins)

    M = generate_spherical_masks(sl);
    num_spheres=length(M);
    sz=size(data_i.tight_liver_mask); 
    [ii,jj,kk] = ind2sub(sz,locations);
    num_maps = length(maps);
    N=(2*sl)+1;

    edges = [linspace(0,1,num_bins+1)];

    for c=1:length(maps)
        patches{c}=zeros(length(ii),(num_spheres*num_bins)); 
    end

    for d=1:length(ii)

        x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,sz(1)),N));
        y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,sz(2)),N));
        z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,sz(3)),N));
        sample_locs = combvec(x,y,z);

        inds = sub2ind(sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));

        for c=1:length(maps)
            for s=1:num_spheres
                h=histcounts(maps{c}(inds(M{s})),edges);
                patches{c}(d,((s-1)*num_bins)+1:s*num_bins) = h;
            end
        end

    end

    temp=[];
    for c=1:num_maps
       temp=[temp,patches{c}];  
    end
end

function s = compute_sensitivity(classification,labels)
    TP = sum(classification.*labels);
    FN = sum(labels.*(1-classification));

    s = TP/(TP+FN);
end

function s = compute_specificity(classification,labels)
    TN = sum((1-classification).*(1-labels));
    FP = sum(classification.*(1-labels));

    s = TN/(TN+FP);
end

function s = compute_accuracy(classification,labels)
    TN = sum((1-classification).*(1-labels));
    FP = sum(classification.*(1-labels));

    s = TN/(TN+FP);
end

function s = compute_precision(classification,labels)
    TP = sum(classification.*labels);
    FP = sum(classification.*(1-labels));

    s = TP/(TP+FP);
end

function DSC = compute_DSC(classification,lab)
    DSC=2*sum(classification.*lab)/(sum(classification)+sum(lab));
end

function M = generate_spherical_masks(sl)
    A=zeros((2*sl)+1,(2*sl)+1,(2*sl)+1);
    center=[sl+1,sl+1,sl+1]; 

    for i=1:size(A,1)
        for j=1:size(A,2)
            for k=1:size(A,3)
                A(i,j,k) = norm([i,j,k] - center);  
            end
        end
    end

    A=round(A); 
    ind=1;
    for r=1:sl
        M{ind} = find(A==r);
        ind=ind+1; 
    end
    return 
end

function maps = generate_maps_dummy(class,sz,locs,num_classes)
    for c=2:num_classes
        maps{c-1,1} = generate_classification_visualization(class(:,c),sz,locs);
    end
end

function map = generate_classification_visualization(p,sz,locs) 
    map=zeros*ones(sz);
    for i=1:length(locs)
        map(locs(i))=p(i); 
    end
    map(1,1,:) = ones(1,size(map,3));
end