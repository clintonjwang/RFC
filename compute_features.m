function [train] = compute_features(data, data_dir, train_indices)

% num_patients = length(data);
% 
% disp('Generating feature cell array...'); 
% tic
% for i=1:length(data)
%     features{i} = generate_feature_array(data{i});
% end
% toc
% 
% disp('Normalizing data...');
% tic
% pr=normalize_data(data,features);
% toc
% 
% data_1 = load(['features/data_',num2str(1),'.mat']);
% data_1 = data_1.data_i;
% 
% num_intensity_maps = length(data_1.p_im);
% num_frangi_maps = length(data_1.frangi); 
% num_grad_maps = length(data_1.grad); 
% num_sf_maps = length(data_1.sf); 
% num_haralick = 3; 
% num_mode_maps = length(data_1.mode);
% 
% clear data_1

tic
for i=1:length(data)
    disp('Generating feature cell array...'); %fast
    features{i} = struct;
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
    features{i} = generate_feature_array(data{i});
end

disp('Normalizing data...');  %fast
[data, pr] = normalize_data(data,features); 
toc

parfor i=1:length(data)
    features{i} = compute_features_single(data{i}, pr);
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


% function: generate_training_data %
function train = generate_training_data(train,features)

load_dir='features';
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
