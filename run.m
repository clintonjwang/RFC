% textscan(fopen('config.txt','r'),'%s','Delimiter','\n')
addpath(genpath('./utils'));
addpath(genpath('./scripts'));
addpath(genpath('./additional'));
train_dir = './raw_imgs';
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

save('./params/train.mat','train','-v7.3');
save('./params/train_indices.mat','train_indices');
save('./params/test_indices.mat','test_indices');

intensities_bin;
%delete(gcp('nocreate'));
tissue_classification;

test_indices = test_indices.test_indices{6};
for idx = test_indices
    disp(idx);
    f = load([feature_dir,'/features_',num2str(idx),'.mat']);
    gt_img = zeros(f.sz);
    pred_img = zeros(f.sz);
    
%     Display predicted image
    for pix_idx = 1:length(f.locations)
        pred_img(f.locations(pix_idx)) = f.labels(pix_idx);
    end
    image(pred_img(:,:,round(f.sz(3)*2/3)), 'CDataMapping','scaled');
    colorbar;
    
%     Display ground truth image
    data = load_nii([train_dir, '/', f.patID, '/nii_files/20s_isotropic.nii.gz']);
    img = double(flip_image(data.img));
    image(img(:,:,round(f.sz(3)*2/3)));
    
%     load([train_data_dir,'/data_',num2str(idx),'.mat']);
end