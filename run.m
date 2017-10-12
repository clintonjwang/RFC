% textscan(fopen('config.txt','r'),'%s','Delimiter','\n')
addpath(genpath('./utils'));
addpath(genpath('./scripts'));
addpath(genpath('./additional'));
train_dir = './training_data/';
train_data_dir = './data';
feature_dir = './features';
param_dir = './params';

patients = dir(train_dir);
patients = {patients.name};
% dir_patients = {dir_patients(4:58).name};
data = acquire_data(patients(3:end), train_dir);

save([train_data_dir,'/data.mat'],'data','-v7.3');
for i = 1:length(data)
    data_i = data{i};
    save([train_data_dir,'/data_',num2str(i),'.mat'],'data_i')   
end

[train,train_indices,test_indices] = compute_features([train_data_dir,'/data.mat'], train_data_dir);

save('./params/train.mat','train');
save('./params/train_indices.mat','train_indices');
save('./params/test_indices.mat','test_indices');

intensities_bin;
tissue_classification;