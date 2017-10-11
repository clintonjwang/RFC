function [train,train_indices,test_indices] = ...
    compute_features(load_data_dir,...
    load_data_pats)

load(load_data_dir); 

num_patients = length(data);

folds=6;
ep=round(linspace(1,num_patients,folds+1));

for i=1:folds
   test_indices{i} = ep(i):ep(i+1); 
   train_indices{i} = setdiff(1:num_patients,test_indices{i}); 
end

disp('Generating feature cell array...'); 
tic
features = generate_feature_array(data,num_patients);
toc

disp('Normalizing data...'); 
tic
pr=normalize_data(data,features); 
toc

data_1 = load([load_data_pats,'/data_',num2str(1),'.mat']);
data_1 = data_1.data_i;

num_intensity_maps = length(data_1.p_im);
num_frangi_maps = length(data_1.frangi); 
num_grad_maps = length(data_1.grad); 
num_sf_maps = length(data_1.sf); 
num_haralick = 3; 
num_mode_maps = length(data_1.mode);

clear data_1

parfor i=1:length(data)
    
    disp(['Current patient... ',num2str(i)]); 
    
    data_i = load([load_data_pats,'/data_',num2str(i),'.mat']);
    data_i = data_i.data_i; 
    
    features{i}.sz = size(data_i.tight_liver_mask); 
    
    disp('Extracting voxel-wise intensities...'); 
    tic
    features{i}.intensities = zeros(length(features{i}.locations),...
        num_intensity_maps);
    
    for j=1:num_intensity_maps
        features{i}.intensities(:,j) = data_i.p_im{j}(features{i}.locations); 
    end
    toc
    
    disp('Extracting mode-map intensities...'); 
    tic
    features{i}.mode_intensities = zeros(length(features{i}.locations),...
        num_mode_maps);
    
    for j=1:num_mode_maps
        features{i}.mode_intensities(:,j) = data_i.mode{j}(features{i}.locations);
    end
    toc
    
    disp('Extracting Frangi intensities...'); 
    tic
    features{i}.frangi = zeros(length(features{i}.locations),...
        num_frangi_maps);
    
    for j=1:num_frangi_maps
        features{i}.frangi(:,j) = data_i.frangi{j}(features{i}.locations); 
    end
    toc
    
    disp('Extracting t2 intensities...');
    tic
    features{i}.t2 = zeros(length(features{i}.locations),1);
    features{i}.t2(:,1) = data_i.t2(features{i}.locations);
    toc
    
    disp('Extracting surface distances...'); 
    tic
    features{i}.intensities = [features{i}.intensities,...
        compute_surface_distance(data_i,features{i}.locations)];
    toc
    
    disp('Extracting gradient intensities...'); 
    tic
    features{i}.gradient = zeros(length(features{i}.locations),...
        num_grad_maps);
    
    for j=1:num_grad_maps
        features{i}.gradient(:,j) = data_i.grad{j}(features{i}.locations);
    end
    toc
    
    disp('Extracting std filter intensities...');
    tic
    features{i}.sf = zeros(length(features{i}.locations),...
        num_sf_maps);
    
    for j=1:num_sf_maps
        features{i}.sf(:,j) = data_i.sf{j}(features{i}.locations);
    end
    toc
    
    disp('Extracting Haralick feature intensities...');
    tic
    features{i}.haralick = zeros(length(features{i}.locations),num_intensity_maps*...
        num_haralick);
    
    ind=1;
    for j1=1:num_intensity_maps
        for j2=1:num_haralick
            features{i}.haralick(:,ind) = data_i.glcm{j1,j2}(...
                features{i}.locations);
            ind=ind+1;
        end
    end
    toc
    
    features{i}.intensities = [features{i}.intensities,...
        features{i}.mode_intensities];
    features{i}=rmfield(features{i},'mode_intensities');
    
    features{i}.intensities = [features{i}.intensities,...
        features{i}.haralick];
    features{i}=rmfield(features{i},'haralick');
    
    features{i}.intensities = [features{i}.intensities,...
        features{i}.sf];
    features{i}=rmfield(features{i},'sf');
    
    features{i}.intensities = [features{i}.intensities,...
        features{i}.gradient];
    features{i}=rmfield(features{i},'gradient');
    
    features{i}.intensities = [features{i}.intensities,...
        features{i}.frangi];
    features{i}=rmfield(features{i},'frangi');
    
    features{i}.intensities = [features{i}.intensities,...
        features{i}.t2];
    features{i}=rmfield(features{i},'t2');
    
    disp('Appending contextual features');
    tic
    features{i} = append_context_features(features{i},data_i,pr);
    toc
    
    features{i}.auto_context_features=[]; 
    features{i}.auto_context_features_boost=[]; 
    
    save_dummy(features{i},i);

    features{i}=features{i}.labels;
  
end 

disp('Extracting training data');  
tic
for model=1:length(train_indices)
    [train.train_labels{model},train.train_patients{model},train.train_locs{model}] = ...
        extract_training_data(features,train_indices{model});
end

train = generate_training_data(train,features);
toc

return
end

function [] =save_dummy(f,i)

%save_dest='C:/Users/John/Desktop/RESULTS/57PATS_July/patient_features';

save_dest='./features';
save([save_dest,'/features_',num2str(i),'.mat'],'f','-v7.3');

return 
end

function features_i = append_context_features(features_i,data_i,pr)

maps{1}=data_i.mode{1}; 
maps{2}=data_i.mode{2};
maps{3}=data_i.mode{3};
maps{4}=data_i.mode{4};

sl=4;

features_i.intensities = [features_i.intensities,...
    compute_patches(maps,features_i,sl,data_i)];

return 
end

function temp = compute_patches(maps,features_i,sl,data_i)
[ii,jj,kk] = ind2sub(features_i.sz,features_i.locations);
num_maps = length(maps);
N=4;

for c=1:length(maps)
    patches{c}=zeros(length(ii),N^3); 
end

for d=1:length(ii)
    
    x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,features_i.sz(1)),N));
    y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,features_i.sz(2)),N));
    z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,features_i.sz(3)),N));
    sample_locs = combvec(x,y,z);
    
    for n=1:(N^3)
        if(data_i.tight_liver_mask(sample_locs(1,n),sample_locs(2,n),...
                sample_locs(3,n))==0)
            sample_locs(:,n)=[ii(d);jj(d);kk(d)];
        end
    end
    
    inds = sub2ind(features_i.sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));

    for c=1:length(maps)
        patches{c}(d,:) = maps{c}(inds); 
    end
    
end

temp=[];
for c=1:num_maps
   temp=[temp,patches{c}];  
end

return
end

% function: generate_training_data %
function train = generate_training_data(train,features)

load_dir='./features';
num_models = length(train.train_labels); 

load([load_dir,'/features_',num2str(1),'.mat']);
nf = size(f.intensities,2);

for model=1:num_models
    num_train_points = length(train.train_patients{model});
    
    train.data{model}=zeros(num_train_points,nf);
    
    for p=1:length(features)
        loc_p = find(train.train_patients{model}==p);
        if(~isempty(loc_p))
            load([load_dir,'/features_',num2str(p),'.mat']);
            for lp=loc_p'
                ind = train.train_locs{model}(lp);
                train.data{model}(lp,:) = f.intensities(ind,:);
            end
        end
    end
end

return
end

% function: generate_feature_array %
function features = generate_feature_array(data,num_patients)

disp('Generating feature struct...');
for i=1:num_patients
    disp(['Patient... ',num2str(i)]);
    features{i} = struct;
    features{i}.patID = data{i}.patID;
    
    features{i}.locations = find(data{i}.tight_liver_mask);
    features{i}.labels=zeros(length(features{i}.locations),1);
    
    for c=1:length(features{i}.locations)
        if(data{i}.necrosis_mask(features{i}.locations(c))==1)
            features{i}.labels(c)=3;
        elseif(data{i}.tumor_mask(features{i}.locations(c))==1)
            features{i}.labels(c)=1;
        elseif(data{i}.vessel_mask(features{i}.locations(c))==1)
            features{i}.labels(c)=2;
        end
    end
end

%{
disp('Generating feature struct...');
parfor i=1:num_patients
    disp(['Patient... ',num2str(i)]);
    features{i} = struct;
    features{i}.patID = data{i}.patID;
    features{i}.locations = find(data{i}.tight_liver_mask);
end
%}
return
end

function temp = compute_surface_distance(data_i,locations)

D=bwdist(1-data_i.tight_liver_mask);

temp=zeros(length(locations),1);
for j=1:length(locations)
    temp(j)=D(locations(j));
end

temp = reshape(temp,[numel(temp),1]); 

return
end

%{
function temp = compute_patches_spherical_distributions(...
    maps,features_i,sl,data_i,pr)

M = generate_spherical_masks(sl);
num_spheres=length(M);
[ii,jj,kk] = ind2sub(features_i.sz,features_i.locations);
num_maps = length(maps);
N=(2*sl)+1;
num_bins=25; 

for c=1:length(maps)
    edges{c} = [-Inf,linspace(pr{c}(1),pr{c}(2),num_bins-1),+Inf];
end

for c=1:length(maps)
    patches{c}=zeros(length(ii),num_bins*num_spheres); 
end

for d=1:length(ii)
    
    x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,features_i.sz(1)),N));
    y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,features_i.sz(2)),N));
    z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,features_i.sz(3)),N));
    sample_locs = combvec(x,y,z);
    
    for n=1:(N^3)
        if(data_i.tight_liver_mask(sample_locs(1,n),sample_locs(2,n),...
                sample_locs(3,n))==0)
            sample_locs(:,n)=[ii(d);jj(d);kk(d)];
        end
    end
    
    inds = sub2ind(features_i.sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));
    
    for c=1:length(maps)
        for s=1:num_spheres
            patches{c}(d,((s-1)*num_bins)+1:s*num_bins) = ...
                histcounts(maps{c}(inds(M{s})),edges{c});
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