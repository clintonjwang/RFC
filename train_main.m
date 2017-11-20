% Add paths, set paths
addpath(genpath('subroutines'));
addpath(genpath('utils'));
addpath(genpath('scripts'));
addpath(genpath('additional'));

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

% Collect images and whole liver masks
patients = dir(train_dir);
patients = {patients.name};
data = acquire_data(patients(3:end), train_dir, train_bool);

save([feature_dir,'/data.mat'],'data','-v7.3');
% for i = 1:length(data)
%     data_i = data{i};
%     save([train_dir,'/data_',num2str(i),'.mat'],'data_i')
% end
% data = load([feature_dir,'/data.mat']); 

[train,train_indices,test_indices] = compute_features(data, train_dir);

save('params/train.mat','train','-v7.3');
save('params/train_indices.mat','train_indices');
save('params/test_indices.mat','test_indices');


train_size = 20;
for i = 1:train_size
    f = intensities_bin(f, i);
    locations{i} = f.locations;
end
save('./params/locations.mat','locations','-v7.3');

f = tissue_classification(train_bool);

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
%     data = load_nii([test_dir, '/', f.patID, '/nii_files/20s_isotropic.nii.gz']);
%     img = double(flip_image(data.img));
%     image(img(:,:,round(f.sz(3)*2/3)));
    
%     load([train_data_dir,'/data_',num2str(idx),'.mat']);
end