function features = compute_features_single(data)
num_patients = length(data);

tic
disp('Generating feature cell array...');
features = generate_feature_array(data,num_patients);
disp('Normalizing data...'); 
pr = normalize_data(data,features); 
toc

data = data.data_i;
num_intensity_maps = length(data.p_im);
num_frangi_maps = length(data.frangi); 
num_grad_maps = length(data.grad); 
num_sf_maps = length(data.sf);
num_mode_maps = length(data.mode);
num_haralick = 3;

clear data_1

% parfor i=1:length(data)

% disp(['Current patient... ',num2str(i)]); 

% data_i = load([load_data_pats,'/data_',num2str(i),'.mat']);

features.sz = size(data.tight_liver_mask); 

disp('Extracting voxel-wise intensities...'); 
tic
features.intensities = zeros(length(features.locations),...
    num_intensity_maps);

for j=1:num_intensity_maps
    features.intensities(:,j) = data.p_im{j}(features.locations); 
end
toc

disp('Extracting mode-map intensities...'); 
tic
features.mode_intensities = zeros(length(features.locations),...
    num_mode_maps);

for j=1:num_mode_maps
    features.mode_intensities(:,j) = data.mode{j}(features.locations);
end
toc

disp('Extracting Frangi intensities...'); 
tic
features.frangi = zeros(length(features.locations),...
    num_frangi_maps);

for j=1:num_frangi_maps
    features.frangi(:,j) = data.frangi{j}(features.locations); 
end
toc

disp('Extracting t2 intensities...');
tic
features.t2 = zeros(length(features.locations),1);
features.t2(:,1) = data.t2(features.locations);
toc

disp('Extracting surface distances...'); 
tic
features.intensities = [features.intensities,...
    compute_surface_distance(data,features.locations)];
toc

disp('Extracting gradient intensities...'); 
tic
features.gradient = zeros(length(features.locations),...
    num_grad_maps);

for j=1:num_grad_maps
    features.gradient(:,j) = data.grad{j}(features.locations);
end
toc

disp('Extracting std filter intensities...');
tic
features.sf = zeros(length(features.locations),...
    num_sf_maps);

for j=1:num_sf_maps
    features.sf(:,j) = data.sf{j}(features.locations);
end
toc

disp('Extracting Haralick feature intensities...');
tic
features.haralick = zeros(length(features.locations),num_intensity_maps*...
    num_haralick);

ind=1;
for j1=1:num_intensity_maps
    for j2=1:num_haralick
        features.haralick(:,ind) = data.glcm{j1,j2}(...
            features.locations);
        ind=ind+1;
    end
end
toc

features.intensities = [features.intensities,...
    features.mode_intensities];
features=rmfield(features,'mode_intensities');

features.intensities = [features.intensities,...
    features.haralick];
features=rmfield(features,'haralick');

features.intensities = [features.intensities,...
    features.sf];
features=rmfield(features,'sf');

features.intensities = [features.intensities,...
    features.gradient];
features=rmfield(features,'gradient');

features.intensities = [features.intensities,...
    features.frangi];
features=rmfield(features,'frangi');

features.intensities = [features.intensities,...
    features.t2];
features=rmfield(features,'t2');

disp('Appending contextual features');
tic
features = append_context_features(features,data,pr);
toc

features.auto_context_features=[]; 
features.auto_context_features_boost=[]; 

save_dummy(features,i);

features=features.labels;




disp('Extracting training data');  
tic
for model=1:length(train_indices)
    [features.train_labels{model},features.train_patients{model},features.train_locs{model}] = ...
        extract_training_data(features,train_indices{model});
end

features = generate_training_data(features,features);
toc

return
end
