% Treilhard J. et al. (2017) Liver Tissue Classification in Patients with
% Hepatocellular Carcinoma by Fusing Structured and Rotationally Invariant
% Context Representation. In: Descoteaux M., Maier-Hein L., Franz A.,
% Jannin P., Collins D., Duchesne S. (eds) Medical Image Computing and
% Computer-Assisted Intervention ? MICCAI 2017. MICCAI 2017. Lecture Notes
% in Computer Science, vol 10435. Springer, Cham
% DOI https://doi.org/10.1007/978-3-319-66179-7_10

% Add paths, set paths
addpath(genpath('subroutines'));
addpath(genpath('utils'));
addpath(genpath('scripts'));
addpath(genpath('additional'));

% uiwait(msgbox(['This program automatically segments a liver into '...
%     'viable tumor, necrosis, vasculature and parenchyma. '...
%     'It requires registered pre-contrast and contrast enhanced abdominal '...
%     'MRIs (in nifti format) along with the whole liver mask (in .ics/'...
%     '.ids format). It outputs binary masks of the 4 classes (also in '...
%     '.ics/.ids format). Only ICS version 1 is supported. '...
%     'The program asks you to select a patient folder to look for the '...
%     'MRIs/masks in. If it cannot find a file automatically, it will '...
%     'prompt you for it.'], 'Random Forest Cascade utility', 'modal'));

% test_dir = uigetdir('', 'Select the folder containing the patient to segment.');
feature_dir = 'features';
param_dir = 'params';
train_bool = false;
save_files = false;

if test_dir == 0
    return
end

% Collect images and whole liver masks
data = acquire_data_single_pat(test_dir, train_bool);
features = compute_features_single([test_dir,'/data.mat'], test_dir, train_bool);

features = intensities_bin(features);
%delete(gcp('nocreate'));
tissue_classification_single(train_bool);

% Display result
% features = load([feature_dir,'/features_',num2str(idx),'.mat']);
pred_img = zeros(features.sz);
for pix_idx = 1:length(features.locations)
    pred_img(features.locations(pix_idx)) = features.labels(pix_idx);
end
image(pred_img(:,:,round(features.sz(3)*2/3)), 'CDataMapping','scaled');
colorbar;