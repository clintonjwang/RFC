function [ status, message, outputVolume ] = LiverSegTrainInterface( inputVolume, inputParams )
%LIVER_SEG_INTERFACE Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT parameters:
%   inputVolume         - Registered pre-contrast, arterial, venous
%   inputParams         - Random forest parameters
% OUTPUT parameters:
%   status              - bool value indicating success or failure
%   message             - Message that may be used to provide information
%                         on failures or success
%   outputMetaDataDict  - Metadata dictionary of the Output Volume
%   outputVolume        - Output Volume

status = false;

inputHeader = inputVolume{1};

%get number of voxels and number of slices from the header
number_of_voxels = header.NumberOfVoxels;

data_dir = inputParams.data_dir;

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
masks = tissue_classification({''}, model_dir, working_dir, working_dir, train_bool, data_dir, out_dir);

if ~save_features
    [~, ~, ~] = rmdir(working_dir, 's');
end

% Display result
mask_names = {'vasculature_mask', 'necrosis_mask', 'viable_tumor_mask'};
mask_display_names = {'vasculature', 'necrosis', 'viable tumor'};




header, metadata, voxel_data;
voxel_data = user_main(false);
header = cell(4,1);
for i = 1:4
    header{i}.HeaderVersion = '20171216';
    header{i}.NumberOfVoxels = size(voxel_data);
    header{i}.VolumeExtent = size(voxel_data);
    header{i}.SliceThickness = VolumeExtent(3);
    header{i}.RescaleSlope = 1;
    header{i}.RescaleIntercept = 0;
    header{i}.TransformationMatrix = zeros(12,1);
    header{i}.Type = 1; %1 for binary, else 0
    header{i}.Description = 'Segmentation';
end

outputVolume{1} = outputHeader;
outputVolume{2} = outputMetaDataDict;
outputVolume{3} = inputVolume{3};

status = true;
message = 'Successfully segmented input';

end

