dir_patients = dir('./training_data');
addpath('./utils/NIfTI_20140122');
% dir_patients = {dir_patients(4:58).name};
data = acquire_data(dir_patients);
save('./Data/data.mat','data','-v7.3');
for i = 1:length(data)
    data_i = data{i};
    save(['./Data/data_',num2str(i),'.mat'],'data_i')   
end
load_data_dir = './Data/data.mat'
load_data_pats = './Data'
[train,train_indices,test_indices] = compute_features(load_data_dir,load_data_pats);
save('./Parameters/train.mat','train');
save('./Parameters/train_indices.mat','train_indices');
save('./Parameters/test_indices.mat','test_indices');

intensities_bin;
tissue_classification;