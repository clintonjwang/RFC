% Add paths, set paths
addpath(genpath('subroutines'));

uiwait(msgbox(['Use this utility to retrain the segmentation tool. '...
    'Be sure you have specified 4 non-overlapping binary masks representing'...
    'viable tumor, necrosis, vasculature, '...
    'and parenchyma segmentations, each in ICS version 1 format.'...
    'Once started, it may take up to several hours per patient.'], 'Random Forest training utility', 'modal'));

train_dir = uigetdir('', 'Select the folder containing all the training data.');
feature_dir = 'features';
param_dir = 'params';
train_bool = true;

if train_dir == 0
    return
end

addpath(genpath('utils'));
addpath(genpath('scripts'));
addpath(genpath('additional'));

% Collect images and whole liver masks
patients = dir(train_dir);
patients = {patients.name};
data = acquire_data(patients(3:end), train_dir, train_bool);

save([feature_dir,'/data.mat'],'data','-v7.3');
% for i = 1:length(data)
%     data_i = data{i};
%     save([train_dir,'/data_',num2str(i),'.mat'],'data_i')
% end
data = load([feature_dir,'/data.mat']); 

[train,train_indices,test_indices] = compute_features(data, train_dir);

save('params/train.mat','train','-v7.3');
save('params/train_indices.mat','train_indices');
save('params/test_indices.mat','test_indices');

intensities_bin;
tissue_classification(train_bool);