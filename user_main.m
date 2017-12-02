% Treilhard J. et al. (2017) Liver Tissue Classification in Patients with
% Hepatocellular Carcinoma by Fusing Structured and Rotationally Invariant
% Context Representation. In: Descoteaux M., Maier-Hein L., Franz A.,
% Jannin P., Collins D., Duchesne S. (eds) Medical Image Computing and
% Computer-Assisted Intervention - MICCAI 2017. MICCAI 2017. Lecture Notes
% in Computer Science, vol 10435. Springer, Cham
% DOI https://doi.org/10.1007/978-3-319-66179-7_10

% Add paths, set paths
addpath(genpath('subroutines'));
addpath(genpath('utils'));
addpath(genpath('scripts'));
addpath(genpath('additional'));

fast_mode = true;
if ~fast_mode
    uiwait(msgbox(['This program automatically segments a liver into '...
        'viable tumor, necrosis, vasculature and parenchyma. '...
        'It requires registered pre-contrast and contrast enhanced abdominal '...
        'MRIs (in nifti format) along with the whole liver mask (in .ics/'...
        '.ids format). It outputs binary masks of the 4 classes (also in '...
        '.ics/.ids format). Only ICS version 1 is supported. '...
        'The program asks you to select a patient folder to look for the '...
        'MRIs/masks in. If it cannot find a file automatically, it will '...
        'prompt you for it.'], 'Random Forest Cascade utility', 'modal'));
end

if fast_mode
    data_dir = '../small_data/0347479';
else
    data_dir = uigetdir('', 'Select the folder containing the patient to segment.');
    if data_dir == 0
        return
    end
end

feature_dir = 'features';
param_dir = 'params';
train_bool = false;
save_files = false;
patients = {data_dir};
%[f,f,f] = fileparts(fn);

model_dir = 'models';
if ~exist(model_dir,'dir')
    model_dir = uigetdir('', 'Select the folder containing all the models.');
    if model_dir == 0
        return
    end
end

working_dir = 'working';
[~, ~, ~] = mkdir(working_dir);

if fast_mode
    out_dir = 'output_masks';
    [~, ~, ~] = mkdir(out_dir);
else
    out_dir = uigetdir('', 'Select a folder to output the binary masks to.');
    if out_dir == 0
        return
    end
end

% Collect images and whole liver masks
data = acquire_data_single_pat(data_dir, train_bool);

if true%exist([working_dir,'/init_features_1.mat'],'file') == 0
    % Initialize labels and locations
    f = struct;
    f.locations = find(data.tight_liver_mask);
    f.labels = zeros(length(f.locations),1);
    save([working_dir,'/init_features_1.mat'], 'f');
    
    % Compute image features for the entire image
    dcell = {data};
    normalize_data(dcell, working_dir);
    clear dcell;
end

% Separate features based on label
% compute_features(patients, working_dir);
data = load([working_dir,'/norm_data_1.mat']);
data = data.data_i;
f = compute_features_single(data, f);
save([working_dir,'/features_1.mat'],'f');

% Save intensities in a separate bin file
if exist([working_dir,'/intensities_1.bin'],'file') == 0
    intensities_bin({'placeholder.....'}, working_dir);
end

% Train random forest model
tissue_classification({'placeholder.....'}, model_dir, working_dir, working_dir, train_bool, data_dir, out_dir);

% Display result
load([working_dir,'/classified_features_2.mat']);
pred_img = zeros(f.sz);
for pix_idx = 1:length(f.locations)
    pred_img(f.locations(pix_idx)) = f.classification{2}(pix_idx);
end
image(pred_img(:,:,round(f.sz(3)*2/3)), 'CDataMapping','scaled');
colorbar;