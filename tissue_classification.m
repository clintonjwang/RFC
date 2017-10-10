function [] = tissue_classification()

dir_main = '/gpfs/home/fas/duncan/jy498/Research/Junlin'
dir_scratch60 = '/gpfs/scratch60/fas/duncan/jy498'
%load relevant global variables 
%load('data_light.mat'); %load relevant global variables 
load([dir_main, '/Data/data.mat']);
load([dir_main, '/Parameters/train_indices.mat']); %holds subset of training patients for each round of cv
load([dir_main, '/Parameters/test_indices.mat']); %holds subset of training patients for each round of cv
load([dir_main, '/Parameters/locations.mat']); %holds subset of training patients for each round of cv
load([dir_main, '/Parameters/train.mat']); %holds subset of training patients for each round of cv

%initialize global variables 
ntrees=1000; %trees in each random forest 
RAC=4; %number of iterations of classification 
sl=8; %structured context patch size 
sl_spherical=5; %spherical context patch size 
.........................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
num_bins=5; %number of histogram bins for spherical histogram context features 
sa=6; %number of auto-context features to sample per dimension 
num_models = length(train.train_labels); %number of folds of cv 
num_patients = length(data); %total number of patients 
min_leaf_size=50; %random forest parameter: mininum leaf size 
op = statset('UseParallel',1); %random forest parameter: train trees in parallel 
%parpool(4); %number of cores to distribute processing over
parpool(12);
num_classes=4; %number of tissue classes 
M = generate_spherical_masks1(sl_spherical); %generate spherical masks 

%compute number of auto-context features...
num_context_features=(sa^3) * (num_classes - 1); %structured context features
num_context_features = num_context_features + ((num_classes-1)*num_bins*length(M)); %plus spherical context features 

disp('Beginning cross-validation...')

%loop through cross-validation
%for model=1:num_models
for model =6
    %loop through rounds of auto-context
    for r=1:RAC
    %for r=4
        
        disp(['Current round of auto-context...',num2str(r)]);
        
        %clear trees from previous round of label-context
        clear C
        %calculate number of predictors (depending on r)
        load([dir_main,'/Feature/features_1.mat']);
        if(r==1)
            num_predictors=f.num_intensity_features;
        else
            num_predictors=f.num_intensity_features+num_context_features;
        end
        clear f
        
        %number of training data points 
        num_data_points = size(train.data{model},1);
        
        %train random forest
        for t=1:ntrees
            %bagging training data 
            bag = randi(num_data_points,[1,num_data_points]);
            if(r==1)
                d=train.data{model}(bag,:);
            else
                d=train_data_ac(bag,:);
            end
            l=train.train_labels{model}(bag,:);
            
            C{t} = compact(fitctree(d,l,'MinLeafSize',min_leaf_size,...
                'NumVariablesToSample',sqrt(num_predictors),...
                'SplitCriterion','deviance','Surrogate','on'));
        end
        % TODO: here to save RF 
        clear train_data_ac;
        
        %compute label probability maps for each tissue class and patient 
        parfor td=1:num_patients
        %for td=1:num_patients
            
            %load corresponding patient 
            f=load([dir_main,'/Feature/features_',num2str(td),'.mat']);
            f=f.f;
            
            disp(['Computing classification for patient...',num2str(td)]);
            
            %read in appearance feature binary file for patient 
            fileID = fopen([dir_scratch60, '/Feature/intensities_',num2str(td),'.bin'],'r');
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
                %TODO: classifier funciton given saved RF
            %for subsequent rounds of classification, incorporate auto-context 
            %features 
            else
                %read in binary file holding auto-context features for
                %patient
                fileID = fopen([dir_scratch60, '/Feature/context_',num2str(td),'.bin'],'r');
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
            if(ismember(td,test_indices{model}))
                
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
                save_dummy(f,td,dir_main);
                
            end
            
        end
        
        %clear the previously trained random forest 
        clear RF
        
        %compute the auto-context features, UNLESS it is the final
        %iteration of classification 
        if(r~=RAC)
            train_data_ac = extract_autocontext_features(...
                classification,data,locations,sl,sa,train.data{model},...
                train.train_patients{model},train.train_locs{model},train_indices{model},...
                num_context_features,sl_spherical,num_bins,dir_main,dir_scratch60);
        end
        
    end
end

return
end






%dummy function to allow variable saves inside parfor-loop 
function [] = save_dummy(f,td,dir_main)

save([dir_main, '/Feature/features_',num2str(td),'.mat'],'f','-v7.3');

return
end

%dummy function to allow binary file saves inside parfor-loop 
function [] = save_dummy_binary(A,td,dir_scratch60)

fileID = fopen([dir_scratch60, '/Feature/context_',num2str(td),'.bin'],'w');
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
    train_locs,train_indices,num_context_features,sl_spherical,num_bins,dir_main,dir_scratch60)

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
%for td=1:num_patients
    
    td
    
    %load patient features
    f=load([dir_main, '/Feature/features_',num2str(td),'.mat']);
    f=f.f;
    
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
    save_dummy_binary(auto_context_features,td,dir_scratch60);
    
    %if patient is "training patient", then save off the relevant auto-context features
    %for future training purposes 
    if(ismember(td,train_indices))
        T{td} = find(train_patients==td);
        acf{td} = auto_context_features(train_locs(T{td}),:);
    end
    
end

%augment the appearance features in the training data with auto-context
%features
temp=zeros(num_train_points,num_context_features);
for td=train_indices
    temp(T{td},:) = acf{td};
end

train_data_ac = [train_data,temp];

return
end



%function to generate 3D probability maps from probability vectors 
function maps = generate_maps_dummy(class,sz,locs,num_classes)

%sigma=[1,2];

for c=2:num_classes
    maps{c-1,1} = generate_classification_visualization(class(:,c),sz,locs);
    %{
    for j=1:numel(sigma)
        maps{c-1,j+1} = imgaussfilt3(maps{c-1,1},sigma(j));
    end
    %}
end
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

%{
function temp = compute_patches_spherical_harmonic(maps,...
    locations,sl,data_i)

[M,Y] = generate_spherical_masks(sl);
num_spheres=size(Y,1);
num_frequencies=size(Y,2);
sz=size(data_i.tight_liver_mask); 
[ii,jj,kk] = ind2sub(sz,locations);
num_maps = length(maps);
N=(2*sl)+1;

for c=1:length(maps)
    patches{c}=zeros(length(ii),(num_spheres*num_frequencies)); 
end

for d=1:length(ii)
    
    x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,sz(1)),N));
    y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,sz(2)),N));
    z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,sz(3)),N));
    sample_locs = combvec(x,y,z);
    
    inds = sub2ind(sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));
    
    for c=1:length(maps)
        for s=1:num_spheres
            local_patch = maps{c}(inds(M{s}));
            for f=1:num_frequencies
                temp = Y{s,f} * local_patch';
                patches{c}(d,((s-1)*num_spheres)+f) = norm(temp);
            end
        end
    end
    
end

temp=[];
for c=1:num_maps
   temp=[temp,patches{c}];  
end

return
end
%}
