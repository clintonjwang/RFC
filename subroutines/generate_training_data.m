function train = generate_training_data(patients, working_dir)
%GENERATE_TRAINING_DATA Summary of this function goes here
%   features is not modified

f = load([working_dir,'/features_',num2str(1),'.mat']);
f = f.f;

% if ~isfield(f, 'num_intensity_features')
nf = size(f.intensities,2);

disp('Extracting training data...');

for i=1:length(patients)
    features{i} = load([working_dir,'/features_',num2str(i),'.mat']);
    features{i} = features{i}.f.labels;
end

train_indices = 1:length(patients);
[train.train_labels,train.train_patients,train.train_locs] = ...
    extract_training_data(features,train_indices);

num_train_points = length(train.train_patients);
train.data=zeros(num_train_points,nf);

for i=1:length(patients)
    loc_i = find(train.train_patients==i);

    if(~isempty(loc_i))
        features{i} = load([working_dir,'/features_',num2str(i),'.mat']);
        features{i} = features{i}.f;

        for li=loc_i'
            ind = train.train_locs(li);
            train.data(li,:) = features{i}.intensities(ind,:);
        end
    end
end

save([working_dir,'/train.mat'], 'train', '-v7.3');
% end

return
end

