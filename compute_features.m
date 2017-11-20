function [train,train_indices,test_indices] = ...
    compute_features(data, data_dir)

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

data_1 = data{1};
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
    
    data_i = load([data_dir,'/data_',num2str(i),'.mat']);
    data_i = data_i.data_i; 
    
    features{i} = compute_features_single(data_i);
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