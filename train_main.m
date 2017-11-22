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

train_bool = true;
model_dir = 'models';
working_dir = 'working_dir';
[~, ~, ~] = mkdir(model_dir);
[~, ~, ~] = mkdir(working_dir);

patients = dir(train_dir);
filenames = {patients.name};
patients = filenames([patients.isdir]);
patients = patients(3:end);
num_patients = length(patients);

% Collect images and whole liver masks
acquire_data(patients, train_dir, working_dir, train_bool);

% Compute features
compute_features(patients, working_dir);

% Generate training data
generate_training_data(patients, working_dir);

% Collect intensities (?)
intensities_bin(patients, working_dir);

% Train random forest model
tissue_classification(patients, model_dir, working_dir, working_dir, train_bool);

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