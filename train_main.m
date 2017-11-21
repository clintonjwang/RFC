% Add paths, set paths
addpath(genpath('subroutines'));
addpath(genpath('utils'));
addpath(genpath('scripts'));
addpath(genpath('additional'));

fast_mode = true;
if ~fast_mode
    uiwait(msgbox(['Use this utility to retrain the segmentation tool. '...
        'All the training data should be stored in a single folder '...
        'where each subfolder contains the data for a single patient. '...
        'There should be no extraneous subfolders. '...
        'Be sure that each patient subfolder has axial precontrast, arterial, ' 
        'venous, and delayed T1-w images, axial T2 images, a binary mask '...
        'representing the whole liver, as well as 4 non-overlapping '...
        'binary masks representing viable tumor, necrosis, vasculature, '...
        'and parenchyma segmentations, each in ICS version 1 format.'...
        'Once started, it may take up to several hours per patient.'], 'Random Forest training utility', 'modal'));
end

if fast_mode
    train_dir = '../small_data';
else
    train_dir = uigetdir('', 'Select the folder containing all the training data.');
    if train_dir == 0
        return
    end
end

feature_dir = 'features';
param_dir = 'params';
train_bool = true;

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

num_patients = length(data);
train_indices = 1:num_patients;
test_indices = [];
train = compute_features(data, train_dir, train_indices);
save('params/train.mat','train','-v7.3');
save('params/train_indices.mat','train_indices');
save('params/test_indices.mat','test_indices');

train_size = length(train_indices);
for i = 1:train_size
    [f{i}, intensities{i}] = intensities_bin(f, i);
    
    fileID = fopen([feat_dir,'/intensities_',num2str(i),'.bin'],'w');
    fwrite(fileID,intensities,'double');
    fclose(fileID);
    save([feat_dir,'/features_',num2str(i),'.mat'],'f','-v7.3'); 
    
    locations{i} = f{i}.locations;
end
save('./params/locations.mat','locations','-v7.3');

f = tissue_classification(train_bool);

% for idx = test_indices
%     disp(idx);
%     f = load([feature_dir,'/features_',num2str(idx),'.mat']);
%     gt_img = zeros(f.sz);
%     pred_img = zeros(f.sz);
%     
% %     Display predicted image
%     for pix_idx = 1:length(f.locations)
%         pred_img(f.locations(pix_idx)) = f.labels(pix_idx);
%     end
%     image(pred_img(:,:,round(f.sz(3)*2/3)), 'CDataMapping','scaled');
%     colorbar;
%     
% %     Display ground truth image
% %     data = load_nii([test_dir, '/', f.patID, '/nii_files/20s_isotropic.nii.gz']);
% %     img = double(flip_image(data.img));
% %     image(img(:,:,round(f.sz(3)*2/3)));
%     
% %     load([train_data_dir,'/data_',num2str(idx),'.mat']);
% end