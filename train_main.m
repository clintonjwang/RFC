function train_main(fast_mode)
%TRAIN_MAIN entry point for training the random forest tissue classifier.

if nargin < 1
    fast_mode = true;
end

addpath(genpath('subroutines'));
addpath(genpath('utils'));
% addpath(genpath('scripts'));
% addpath(genpath('additional'));

if ~fast_mode
    uiwait(msgbox(['All training data should be stored in a single folder '...
        'where each subfolder contains the data for a single patient. '...
        'There should be no extraneous subfolders. '...
        'Be sure that each patient subfolder has axial precontrast, arterial, ' 
        'venous, and delayed T1-w images, axial T2 images, a binary mask '...
        'representing the whole liver, as well as 4 non-overlapping '...
        'binary masks representing viable tumor, necrosis, vasculature, '...
        'and parenchyma segmentations, each in ICS version 1 format.'...
        'Once started, it may take over an hour per patient.'], 'Random Forest training', 'modal'));
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
working_dir = 'working_train';
[~, ~, ~] = mkdir(model_dir);
[~, ~, ~] = mkdir(working_dir);

patients = dir(train_dir);
filenames = {patients.name};
patients = filenames([patients.isdir]);
patients = patients(3:end);
num_patients = length(patients);

% Collect images and whole liver masks
acquire_data(patients, train_dir, working_dir, train_bool);

if exist([working_dir,'/init_features_',num2str(1),'.mat'],'file') == 0
    % Initialize labels and locations
    [data] = init_features(patients, working_dir);
    
    % Compute image features for the entire image
    normalize_data(data, working_dir);
end

% Separate features based on label
compute_features(patients, working_dir);

% Generate training data
if exist([working_dir,'/train.mat'],'file') == 0
    generate_training_data(patients, working_dir);
end

% Save intensities in a separate bin file
if exist([working_dir,'/intensities_1.bin'],'file') == 0    
    intensities_bin(patients, working_dir);
end

% Train random forest model
tissue_classification(patients, model_dir, working_dir, working_dir, train_bool, train_dir);

[~, ~, ~] = rmdir(working_dir, 's');
if ~fast_mode
    uiwait(msgbox(['Training complete. The tree weights have been saved to ',...
        model_dir, '.'], 'Random Forest training complete', 'modal'));
end

end