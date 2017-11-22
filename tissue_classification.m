function tissue_classification(patients, model_dir, working_dir, scratch_dir, train_mode)
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

%loop through rounds of auto-context
for r=1:RAC
    disp(['Current round of auto-context...',num2str(r)]);

    %clear trees from previous round of label-context
    clear C
    
    %calculate number of predictors (depending on r)
    f = load([working_dir,'/features_1.mat']);
    f = f.f;
    if(r==1)
        num_predictors=f.num_intensity_features;
    else
        num_predictors=f.num_intensity_features+num_context_features;
    end
    clear f

    train = load([working_dir,'/train.mat']);
    train = train.train;
    if train_mode
        %number of training data points 
        num_data_points = size(train.data,1);

        %train random forest
        for t=1:ntrees
            %bagging training data 
            bag = randi(num_data_points,[1,num_data_points]);
            if(r==1)
                data=train.data(bag,:);
            else
                data=train_data_ac(bag,:);
            end
            labs=train.train_labels(bag,:);

            C{t} = compact(fitctree(data,labs,'MinLeafSize',min_leaf_size,...
                'NumVariablesToSample',round(sqrt(num_predictors)),...
                'SplitCriterion','deviance','Surrogate','on')); %opt
        end
        
        save([model_dir,'/tree_',num2str(r),'.mat'],'C');
        clear train_data_ac;
        
    else
%         fileID = fopen([model_dir,'/tree_',num2str(r),'.mat'],'r');
%         C = fread(fileID,[numel(f.labels), num_context_features],'double');
%         fclose(fileID);
        load([model_dir,'/tree_',num2str(r),'.mat']);
    end

    %compute label probability maps for each tissue class and patient 
    for td=1:num_patients %parfor
        %load corresponding patient
        f=load([working_dir,'/features_',num2str(td),'.mat']);
        f=f.f;

        disp(['Computing classification for patient...',num2str(td)]);

        %read in appearance feature binary file for patient 
        fileID = fopen([scratch_dir, '/intensities_',num2str(td),'.bin'],'r');
        intensities = fread(fileID,[numel(f.labels),...
            f.num_intensity_features],'double');
        fclose(fileID);

        %variable holding tissue probability estimates 
        classification{td}=zeros(numel(f.labels),num_classes);

        tic
        %if first round of classification: appearance features only 
        if(r==1)
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
        toc

        %if patient is a "test patient" for this round of
        %cross-validation, then compute accuracy measures of the 
        %classification 
        if ~train_mode
            f.classification{r}=classification{td};

            if(r==1)
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
            r
            save_dummy(f,td,working_dir);
        end
    end

    %clear the previously trained random forest 
    clear RF

    %compute the auto-context features, unless it is the final iteration
    if(r~=RAC)
        train_data_ac = extract_autocontext_features(...
            classification,data,locations,sl,sa,train.data,...
            train.train_patients,train.train_locs,train_indices,...
            num_context_features,sl_spherical,num_bins,scratch_dir);
        
    else
        [vasc_mask, nec_mask, viatumor_mask, paren_mask] = get_masks(classification)
        write_ids_mask(vasc_mask, '.', '.', 'vasculature_mask');
        write_ids_mask(nec_mask, '.', '.', 'necrosis_mask');
        write_ids_mask(viatumor_mask, '.', '.', 'viable_tumor_mask');
        write_ids_mask(paren_mask, '.', '.', 'parenchyma_mask');
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

%function to compute accuracy measures of the classification 
function [AUC,sensitivity,specificity,precision,accuracy,DSC]...
    = compute_effectivness(scores,labels)

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
    train_locs,num_context_features,sl_spherical,num_bins,save_dir)

%initialize variables
num_patients = length(classification);
num_classes = 4;
num_train_points = length(train_patients);

disp('Extracting auto-context features...')

for td=1:num_patients
    T{td}=[];
    acf{td}=[];
end

%loop through patients
parfor td=1:num_patients
    disp(td);
    
    %load patient features
    f=load(['features/features_',num2str(td),'.mat']);
    
    %compute label context features
    tic
    maps = generate_maps_dummy(classification{td},f.sz,locations{td},num_classes);
    
    auto_context_features = [compute_patches(maps,locations{td},sl,data{td},sa),...
        compute_patches_spherical(maps,locations{td},sl_spherical,data{td},num_bins)];
    
    %{
    auto_context_features = [compute_patches_spherical(maps,locations{td},...
        sl_spherical,data{td},num_bins)];
    %}
    toc
    
    %save auto-context features in binary file 
    save_dummy_binary(auto_context_features,td,save_dir);
    
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
