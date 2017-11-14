% Add paths, set paths
uiwait(msgbox(['Use this utility to retrain the segmentation tool. '...
    'Be sure you have specified 4 non-overlapping binary masks representing'...
    'viable tumor, necrosis, vasculature, '...
    'and parenchyma segmentations, each in ICS version 1 format.'...
    'Once started, it will take several hours per patient.'], 'qEASLy utility', 'modal'));

train_data_dir = uigetdir('', 'Select the folder containing all the training data.');

if train_data_dir == 0
    return
end

addpath(genpath('utils'));
addpath(genpath('scripts'));
addpath(genpath('additional'));

train_dir = 'raw_imgs/';
feature_dir = 'features';
param_dir = 'params';

% Collect images and whole liver masks
patients = dir(train_dir);
patients = {patients.name};
data = acquire_data(patients(3:end), train_dir);

save([train_data_dir,'/data.mat'],'data','-v7.3');
for i = 1:length(data)
    data_i = data{i};
    save([train_data_dir,'/data_',num2str(i),'.mat'],'data_i')   
end

[train,train_indices,test_indices] = compute_features([train_data_dir,'/data.mat'], train_data_dir);

save('params/train.mat','train','-v7.3');
save('params/train_indices.mat','train_indices');
save('params/test_indices.mat','test_indices');

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