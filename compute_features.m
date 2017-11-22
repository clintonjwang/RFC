function [train, features] = compute_features(data, data_dir, train_indices)

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

compute_features_single('data_1.mat', pr, features{1});
% parfor i=1:length(data)
%     features{i} = compute_features_single('data_1.mat', pr);
% end

disp('Extracting training data');  
tic
[train.train_labels,train.train_patients,train.train_locs] = ...
    extract_training_data(features,train_indices);
save('params/features.mat','features');

train = generate_training_data(train, './features');
toc

return
end
