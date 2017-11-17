% Treilhard J. et al. (2017) Liver Tissue Classification in Patients with
% Hepatocellular Carcinoma by Fusing Structured and Rotationally Invariant
% Context Representation. In: Descoteaux M., Maier-Hein L., Franz A.,
% Jannin P., Collins D., Duchesne S. (eds) Medical Image Computing and
% Computer-Assisted Intervention ? MICCAI 2017. MICCAI 2017. Lecture Notes
% in Computer Science, vol 10435. Springer, Cham
% DOI https://doi.org/10.1007/978-3-319-66179-7_10

% Add paths, set paths
addpath(genpath('subroutines'));

uiwait(msgbox(['This program automatically segments a liver into '...
    'viable tumor, necrosis, vasculature and parenchyma. '...
    'It requires registered pre-contrast and contrast enhanced abdominal '...
    'MRIs (in nifti format) along with the whole liver mask (in .ics/'...
    '.ids format). It outputs binary masks of the 4 classes (also in '...
    '.ics/.ids format). Only ICS version 1 is supported. '...
    'The program asks you to select a patient folder to look for the '...
    'MRIs/masks in. If it cannot find a file automatically, it will '...
    'prompt you for it.'], 'Random Forest Cascade utility', 'modal'));

test_dir = uigetdir('', 'Select the folder containing the patient to segment.');
feature_dir = 'features';
param_dir = 'params';
train_bool = false;
save_files = false;

if test_dir == 0
    return
end

addpath(genpath('utils'));
addpath(genpath('scripts'));
addpath(genpath('additional'));

% Collect images and whole liver masks
data = acquire_data_single_pat(test_dir, train_bool);
features = compute_features_single([test_dir,'/data.mat'], test_dir, train_bool);

intensities_bin;
%delete(gcp('nocreate'));
tissue_classification(train_bool);

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