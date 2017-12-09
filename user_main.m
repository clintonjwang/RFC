function user_main(fast_mode)
%USER_MAIN entry point for using the trained random forest classifier

if nargin < 1
    fast_mode = true;
end

addpath(genpath('subroutines'));

if ~fast_mode
    uiwait(msgbox(['Using this model requires that it has been previously '...
        'trained. The weights should be saved as "tree_x.mat". '...
        'The model also requires registered contrast enhanced T1/T2 '...
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

use_bias_field = true;
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

working_dir = 'working_test';
[~, ~, ~] = mkdir(working_dir);

if fast_mode
    out_dir = [data_dir,'/output_masks'];
    [~, ~, ~] = mkdir(out_dir);
else
    out_dir = uigetdir('', 'Select a folder to output the binary masks to.');
    if out_dir == 0
        return
    end
end

if ~fast_mode
    prompt = {'Save features at the end of the run?',...
            'Do images have separate T1-w bias field corrections saved as a nifti?'};
    dlg_title = 'Run options';
    num_lines = 1;
    defaultans = {'no','yes'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    if answer == 0
        return
    end
else
    answer = {'yes','yes'};
end

save_features = strcmp(answer{1},'yes') == 1;
use_bias_field = strcmp(answer{2},'yes') == 1;


% Collect images and whole liver masks
data = acquire_data_single_pat(data_dir, train_bool, use_bias_field);

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
    intensities_bin({''}, working_dir);
end

% Train random forest model
tissue_classification({''}, model_dir, working_dir, working_dir, train_bool, data_dir, out_dir);

if ~save_features
    [~, ~, ~] = rmdir(working_dir, 's');
end

% Display result
mask_names = {'vasculature_mask', 'necrosis_mask', 'viable_tumor_mask'};
mask_display_names = {'vasculature', 'necrosis', 'viable tumor'};
display_scrolling_mask('20s', data_dir, out_dir, mask_names, mask_display_names);
end