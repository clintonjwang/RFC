function tissue_classification(patients, model_dir, working_dir, scratch_dir, train_bool, data_dir, out_dir)
%TISSUE_CLASSIFICATION Multiple
%   Detailed explanation goes here

%initialize global variables
ntrees=800; %trees in each random forest 
RAC=2; %3%number of iterations of classification
sl=8; %structured context patch size
sl_spherical=5; %spherical context patch size
num_bins=5; %number of histogram bins for spherical histogram context features 
sa=6; %number of auto-context features to sample per dimension
% num_cv_iterations = length(train.train_labels); %number of folds of cv 
num_patients = length(patients); %total number of patients 
min_leaf_size=50; %random forest parameter: mininum leaf size 
% opt = statset('UseParallel',4); %random forest parameter: train trees in parallel 
% parpool(4); %number of cores to distribute processing over
num_classes=4; %number of tissue classes
M = generate_spherical_masks(sl_spherical); %generate spherical masks
num_context_features = (sa^3) * (num_classes - 1)... %structured context features
         + ((num_classes-1)*num_bins*length(M)); %plus spherical context features 

for i=1:num_patients
    data{i} = load([working_dir,'/norm_data_',num2str(i),'.mat']);
    data{i} = data{i}.data_i;
end
    
%loop through rounds of auto-context
for r=1:RAC
    %clear trees from previous round of label-context
    clear C
    
    %calculate number of predictors (depending on r)
    f = load([working_dir,'/small_features_1.mat']);
    f = f.f;
    if(r==1)
        num_predictors=f.num_intensity_features;
    else
        num_predictors=f.num_intensity_features+num_context_features;
    end
    clear f

    if train_bool
        disp(['Fitting random forest (round ', num2str(r), ' of auto-context)...']);

        tic
        train = load([working_dir,'/train.mat']);
        train = train.train;
        %number of training data points 
        num_data_points = size(train.data,1);

        %train random forest
        for t=1:ntrees
            %bagging training data 
            bag = randi(num_data_points,[1,num_data_points]);
            if(r==1)
                tree_data=train.data(bag,:);
            else
                tree_data=train_data_ac(bag,:);
            end
            labs=train.train_labels(bag,:);
            
            C{t} = compact(fitctree(tree_data,labs,'MinLeafSize',min_leaf_size,...
                'NumVariablesToSample',sqrt(num_predictors),...
                'SplitCriterion','deviance','Surrogate','on')); %opt
        end
        
        save([model_dir,'/tree_',num2str(r),'.mat'],'C');
        clear train_data_ac;
        toc
        
    else
%         fileID = fopen([model_dir,'/tree_',num2str(r),'.mat'],'r');
%         C = fread(fileID,[numel(f.labels), num_context_features],'double');
%         fclose(fileID);
        load([model_dir,'/tree_',num2str(r),'.mat']);
    end
    
    disp('Computing classification...');
    
    tic
    %compute label probability maps for each tissue class and patient 
    for td=1:num_patients %parfor
        %load corresponding patient
        f=load([working_dir,'/small_features_',num2str(td),'.mat']);
        f=f.f;
        
        if train_bool
            disp(['Computing classification for patient...',num2str(td)]);
        end
        
        %read in appearance feature binary file for patient 
        fileID = fopen([scratch_dir, '/intensities_',num2str(td),'.bin'],'r');
        intensities = fread(fileID,[numel(f.labels),...
            f.num_intensity_features],'double');
        fclose(fileID);

        %variable holding tissue probability estimates 
        classification{td}=zeros(numel(f.labels),num_classes);

        %if first round of classification: appearance features only 
        if r==1
            for t=1:ntrees
                [~,scores] = predict(C{t},intensities);
                classification{td} = classification{td} + scores;
            end
        %for subsequent rounds of classification, incorporate auto-context 
        %features 
        else
            %read in binary file holding auto-context features for
            %patient
            fileID = fopen([scratch_dir, '/context_',num2str(td),'.bin'],'r');
            auto_context_features = fread(fileID,[numel(f.labels),...
                num_context_features],'double');
            fclose(fileID);
            for t=1:ntrees
                [~,scores] = predict(C{t},[intensities,auto_context_features]);
                classification{td} = classification{td} + scores;
            end
        end
        classification{td}=classification{td}./ntrees;
        
        %if patient is a "test patient" for this round of
        %cross-validation, then compute accuracy measures of the 
        %classification
        if ~train_bool
            f.classification{r}=classification{td};

            if r==1
                f.AUC=[];
                f.sensitivity = [];
                f.specificity = [];
                f.precision = [];
                f.accuracy = [];
                f.DSC = [];
            end

            [AUC,sensitivity,specificity,precision,accuracy,DSC] = ...
                compute_effectivness(f.classification{r},f.labels);
            f.AUC(:,r) = AUC;
            f.sensitivity(:,r) = sensitivity;
            f.specificity(:,r) = specificity;
            f.precision(:,r) = precision;
            f.accuracy(:,r) = accuracy;
            f.DSC(:,r) = DSC;
            %save_dummy(f,td,working_dir);
            save([working_dir,'/classified_features_',num2str(r),'.mat'],'f','-v7.3');
        end
    end
    toc

    %clear the previously trained random forest 
    clear RF
    load([working_dir,'/locations.mat']);
    
    %compute the auto-context features, unless it is the final iteration
    train_indices = 1:num_patients;
    if train_bool
        if r~=RAC
            train_data_ac = extract_autocontext_features(...
                classification,data,locations,sl,sa,train.data,...
                train.train_patients,train.train_locs,...
                num_context_features,sl_spherical,num_bins,scratch_dir);
        end
    else
        if r~=RAC
            extract_autocontext_features_test(...
                classification,data,locations,sl,sa,...
                sl_spherical,num_bins,scratch_dir);
        else
            for td=1:num_patients %parfor
                load([working_dir,'/small_features_',num2str(td),'.mat']);

                scores = classification{td};
                pred_labels = (scores(:,2)>scores(:,1)) .* (scores(:,2)>scores(:,3)).* (scores(:,2)>scores(:,4)) + ...
                     2 * (scores(:,3)>scores(:,1)) .* (scores(:,3)>scores(:,2)).* (scores(:,3)>scores(:,4)) + ...
                     3 * (scores(:,4)>scores(:,1)) .* (scores(:,4)>scores(:,2)).* (scores(:,4)>scores(:,3));
                pred_labels = pred_labels + 1;

                pred_img = zeros(f.sz);
                for pix_idx = 1:length(f.locations)
                    pred_img(f.locations(pix_idx)) = pred_labels(pix_idx);
                end
                rescaled_img = round(imresize3(pred_img, [260 320 88]));
%                 rescaled_img = pred_img;
                vasc_mask = rescaled_img == 2;
                nec_mask = rescaled_img == 3;
                viatumor_mask = rescaled_img == 4;
                paren_mask = rescaled_img == 1;

                [~,~,~]=mkdir(out_dir);
    %             [vasc_mask, nec_mask, viatumor_mask, paren_mask] = get_masks(classification{td});
                write_ids_mask(vasc_mask, data_dir, out_dir, 'vasculature_mask');
                write_ids_mask(nec_mask, data_dir, out_dir, 'necrosis_mask');
                write_ids_mask(viatumor_mask, data_dir, out_dir, 'viable_tumor_mask');
                write_ids_mask(paren_mask, data_dir, out_dir, 'parenchyma_mask');
            end
        end
    end
end

return
end




%dummy function to allow variable saves inside parfor-loop 
function [] = save_dummy(f,td,dir_main)

save(['./features/features_',num2str(td),'.mat'],'f','-v7.3');

return
end

%dummy function to allow binary file saves inside parfor-loop 
function [] = save_dummy_binary(A,td,dir_scratch60)

fileID = fopen([dir_scratch60, '/context_',num2str(td),'.bin'],'w');
fwrite(fileID,A,'double');
fclose(fileID);

return
end

function [AUC,sensitivity,specificity,precision,accuracy,DSC]...
    = compute_effectivness(scores,labels)
% compute accuracy measures of the classification 

%loop through labels
for c=1:3
    %check if patient has corresponding tissue class
    if(length(find(labels==c))>0)
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

return
end



%function to compute auto-context features 
function train_data_ac = extract_autocontext_features(classification,...
    data,locations,sl,sa,train_data,train_patients,...
    train_locs,num_context_features,sl_spherical,num_bins,working_dir)

disp('Extracting auto-context features...')

tic
%initialize variables
num_patients = length(classification);
num_classes = 4;
num_train_points = length(train_patients);

for td=1:num_patients
    T{td}=[];
    acf{td}=[];
end

%loop through patients
for td=1:num_patients %parfor
    %load patient features
    load([working_dir,'/small_features_',num2str(td),'.mat']);
    
    %compute label context features
    maps = generate_maps_dummy(classification{td},f.sz,locations{td},num_classes);
    
    auto_context_features = [compute_patches(maps,locations{td},sl,data{td},sa),...
        compute_patches_spherical(maps,locations{td},sl_spherical,data{td},num_bins)];
    
    %save auto-context features in binary file 
    save_dummy_binary(auto_context_features,td,working_dir);
    
    %if patient is "training patient", then save the relevant auto-context features
    %for future training purposes 
    T{td} = find(train_patients==td);
    acf{td} = auto_context_features(train_locs(T{td}),:);
end

%augment the appearance features in the training data with auto-context
%features
temp=zeros(num_train_points,num_context_features);
for td=1:num_patients
    temp(T{td},:) = acf{td};
end

train_data_ac = [train_data,temp];
toc

return
end

%function to compute auto-context features 
function extract_autocontext_features_test(classification,...
    data,locations,sl,sa,sl_spherical,num_bins,working_dir)

disp('Extracting auto-context features...')
tic

%initialize variables
num_classes = 4;

td=1;
load([working_dir,'/small_features_1.mat']);

%compute label context features
maps = generate_maps_dummy(classification{td},f.sz,locations{td},num_classes);

auto_context_features = [compute_patches(maps,locations{td},sl,data{td},sa),...
    compute_patches_spherical(maps,locations{td},sl_spherical,data{td},num_bins)];

%save auto-context features in binary file 
%     save_dummy_binary(auto_context_features,td,working_dir);
fileID = fopen([working_dir, '/context_1.bin'],'w');
fwrite(fileID, auto_context_features, 'double');
fclose(fileID);

toc

return
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

M = generate_spherical_masks1(sl);
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

return
end

function s = compute_sensitivity(classification,labels)

TP = sum(classification.*labels);
FN = sum(labels.*(1-classification));

s = TP/(TP+FN);

return
end

function s = compute_specificity(classification,labels)

TN = sum((1-classification).*(1-labels));
FP = sum(classification.*(1-labels));

s = TN/(TN+FP);

return
end

function s = compute_accuracy(classification,labels)

TN = sum((1-classification).*(1-labels));
FP = sum(classification.*(1-labels));

s = TN/(TN+FP);

return
end

function s = compute_precision(classification,labels)

TP = sum(classification.*labels);
FP = sum(classification.*(1-labels));

s = TP/(TP+FP);

return
end

function DSC = compute_DSC(classification,lab)

DSC=2*sum(classification.*lab)/(sum(classification)+sum(lab));

return
end
