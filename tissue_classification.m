function masks = tissue_classification(patients, model_dir, working_dir, train_bool, data_dir, mask_dir, config)
%TISSUE_CLASSIFICATION
%   Detailed explanation goes here

    if train_bool
        disp('Training random forest...');
        save([model_dir, '/tree_config.mat'], 'config');
    else
        disp('Traversing random forest...');
        load([model_dir, '/tree_config.mat']);
    end

    ntrees = config.ntrees; %trees in each random forest, 800 by default
    RAC = config.RAC; %number of iterations of classification, 3 by default
    sl=config.sl; %structured context patch size, 8 by default
    sl_spherical=config.sl_spherical; %spherical context patch size, 5 by default
    num_bins=config.num_bins; %number of histogram bins for spherical histogram context features, 5 by default
    sa=config.sa; %number of auto-context features to sample per dimension, 6 by default
    min_leaf_size = config.min_leaf_size; %random forest parameter: mininum leaf size, 50 by default

    num_patients = length(patients); %total number of patients
    num_classes = 4; %number of tissue classes

    %initialize global variables
    opt = statset('UseParallel',true); %random forest parameter: train trees in parallel 
%     parpool(4); %number of cores to distribute processing over
    M = generate_spherical_masks(sl_spherical); %generate spherical masks
    num_context_features = (sa^3) * (num_classes - 1)... %structured context features
             + ((num_classes-1)*num_bins*length(M)); %plus spherical context features 

    train = load_wrapper([working_dir,'/train.mat']);
    
    start_r = 2;
    for td=1:num_patients
        if ~exist([working_dir, '/context_',patients{td},'_1.bin'], 'file')
            start_r = 1;
            break
        end
    end
    
    %loop through rounds of auto-context
    for r=start_r:RAC
        %clear trees from previous round of label-context
        clear C

        if train_bool && ~exist([model_dir,'/tree_',num2str(r),'.mat'], 'file')
            disp(['Fitting random forest (round ', num2str(r), ' of auto-context)...']);

            %calculate number of predictors (depending on r)
            tmp = load_wrapper([working_dir,'/small_features_',patients{1},'.mat']);
            if(r==1)
                num_predictors=tmp.num_intensity_features;
            else
                num_predictors=tmp.num_intensity_features+num_context_features;
            end
            clear tmp
            
            tic
            %train random forest
            for t=1:ntrees
                %bagging training data 
                num_data_points = round(size(train.data,1) / 3);
                bag = randi(size(train.data,1),[1,num_data_points]);
                if(r==1)
                    tree_data=train.data(bag,:);
                else
                    indices = unique(bag);
                    for i=1:length(indices)
                        train_data_subset{i} = train.data(indices(i),:);
                        acf{i} = load_context(patients{indices(i)}, r, num_context_features, working_dir);
                    end
                    if ~exist('train_data_ac', 'var')
                        parfor td=1:num_patients
                            acf{td} = load_context(patients{td}, r, num_context_features, working_dir);
                        end
                        tree_data = [train_data_subset, acf];
                    end
                    tree_data=train_data_ac(bag,:);
                end
                labs=train.train_labels(bag,:);

                C{t} = compact(fitctree(tree_data,labs,'MinLeafSize',min_leaf_size,...
                    'NumVariablesToSample',sqrt(num_predictors),...
                    'SplitCriterion','deviance','Surrogate','on')); %opt
            end

            save_wrapper(C, [model_dir,'/tree_',num2str(r),'.mat']);
            clear train_data_ac;
            toc
        else
            C = load_wrapper([model_dir,'/tree_',num2str(r),'.mat']);
        end

        disp('Computing classification...');
        tic
        %compute label probability maps for each tissue class and patient 
        parfor td=1:num_patients
            %load corresponding patient
            f = load_wrapper([working_dir,'/small_features_',patients{td},'.mat']);

            %read in appearance feature binary file for patient 
            fileID = fopen([working_dir, '/intensities_',patients{td},'.bin'],'r');
            intensities = fread(fileID,[numel(f.labels),...
                f.num_intensity_features],'double');
            fclose(fileID);

            %initialize
            classification{td}=zeros(numel(f.labels),num_classes);

            %if first round of classification: appearance features only 
            if r==1
                for t=1:ntrees
                    [~,scores] = predict(C{t},intensities);
                    classification{td} = classification{td} + scores;
                end
            %for subsequent rounds, incorporate auto-context features 
            else
                fileID = fopen([working_dir, '/context_',patients{td},'_',num2str(r-1),'.bin'],'r');
                auto_context_f = fread(fileID,[numel(f.labels),...
                    num_context_features],'double');
                fclose(fileID);
                for t=1:ntrees
                    [~,scores] = predict(C{t},[intensities,auto_context_f]);
                    classification{td} = classification{td} + scores;
                end
            end
            classification{td}=classification{td}./ntrees;
        end
        toc
        clear RF

        %compute the auto-context features, unless it is the final iteration
        if r<RAC
            if train_bool
                train_data_ac = extract_autocontext_features(...
                    classification,sl,sa, patients,...
                    train.train_locs, num_context_features,...
                    sl_spherical,num_bins,working_dir,r, train_bool);
                train_data_ac = [train.data, train_data_ac];
            else
                extract_autocontext_features(...
                    classification,sl,sa,patients,[],[],...
                    sl_spherical,num_bins,working_dir,r, train_bool);
            end
        end
    end
    
    %% produce masks
    if ~train_bool
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

%function to compute auto-context features 
function acf = extract_autocontext_features(classification,...
    sl, sa, patients, train_locs, num_context_features,...
    sl_spherical,num_bins,working_dir, r, train_bool)

    disp('Extracting auto-context features...')

    tic
    %initialize variables
    num_patients = length(classification);
    num_classes = 4;
    locations = load_wrapper([working_dir,'/locations.mat']);
    
    %loop through patients
    acf = cell(1,num_patients);
    parfor td=1:num_patients
        if exist([working_dir, '/context_',patients{td},'_',num2str(r),'.bin'], 'file')
            acf{td} = load_context(patients{td}, r, num_context_features, working_dir);
            continue
        end
            
        %load patient features
        f = load_wrapper([working_dir,'/small_features_',patients{td},'.mat']);
        data = load_wrapper([working_dir,'/norm_data_',patients{td},'.mat']);
        
        %compute label context features
        maps = generate_maps_dummy(classification{td},f.sz,locations{td},num_classes);

        acf{td} = [compute_patches(maps,locations{td},sl,data,sa),...
            compute_patches_spherical(maps,locations{td},sl_spherical,data,num_bins)];

        %save auto-context features in binary file 
        save_dummy_binary(acf{td}, [working_dir, '/context_',patients{td},'_',num2str(r),'.bin']);
    end
    toc
end

function acf = load_context(patient, r, num_context_features, working_dir)
    f = load_wrapper([working_dir,'/small_features_',patient,'.mat']);
    fileID = fopen([working_dir, '/context_',patient,'_',num2str(r),'.bin'],'r');
    acf = fread(fileID,[numel(f.labels), num_context_features],'double');
    fclose(fileID);
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