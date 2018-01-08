function [ status, message, outputVolumes ] = LiverSegUsageInterface( inputVolumes, inputParams )
%LIVER_SEG_INTERFACE Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT parameters:
%   inputVolume         - Registered pre-contrast, arterial, venous, T1 bias
%   and T2 in that order.
%   inputParams         - Random forest parameters
% OUTPUT parameters:
%   status              - bool value indicating success or failure
%   message             - Message that may be used to provide information
%                         on failures or success
%   outputVolumes       - Output binary masks (viable tumor, necrosis,
%                         vasculature, and parenchyma in that order)

status = false;
train_bool = false;

data = struct;
data.pre = inputVolumes{1}{3};
data.art = inputVolumes{2}{3};
data.pv = inputVolumes{3}{3};
data.t2 = inputVolumes{4}{3};
data.bf = inputVolumes{5}{3};
data.liver_mask = inputVolumes{6}{3};
data.orig_dims = size(data.art);
data.iso_full_size = data.orig_dims;
% make it isotropic later
% make_isotropic_niis(data_dir, R, train_bool, use_bias_field);
% data = load_niis(data, data_dir, train_bool, use_bias_field);

%get T1 image dimensions
[N1,N2,N3] = size(data.art);

%shrink the masks to proper size if necessary
data = shrink_masks(data, N1, N2, N3, train_bool);

%compute tightened liver mask
r=4;
rs=r/data.resolution(1);
D=bwdist(1-data.liver_mask);
data.tight_liver_mask = zeros(N1,N2,N3);
data.tight_liver_mask(D>rs) = data.liver_mask(D>rs);

%find edges of the liver and crop the images to this dimension
[i,j,k] = ind2sub(size(data.liver_mask),find(data.liver_mask));
i_min = min(i);
i_max = max(i);

j_min=min(j);
j_max=max(j);

k_min=min(k);
k_max=max(k);

data.cutoffs_high = [i_max, j_max, k_max];
data.cutoffs_low = [i_min, j_min, k_min];

sf={'pre','art','pv','t2','liver_mask','tumor_mask','bf',...
    'tight_liver_mask','vessel_mask','necrosis_mask'};

for sf_count=1:length(sf)
    if(isfield(data,sf{sf_count}))
        data.(sf{sf_count}) = data.(sf{sf_count})(i_min:i_max, ...
            j_min:j_max,k_min:k_max);
    end
end

%get number of voxels and number of slices from the header
data = inputVolumes;
model_dir = inputParams.model_dir;
working_dir = inputParams.working_dir;
% out_dir = inputParams.out_dir;
save_features = inputParams.save_features;
use_bias_field = isfield(inputVolumes,'T1_bias');

% Initialize labels and locations
f = struct;
f.locations = find(data.tight_liver_mask);
f.labels = zeros(length(f.locations),1);
save([working_dir,'/init_features_1.mat'], 'f');

% Compute image features for the entire image
dcell = {data};
normalize_data(dcell, working_dir);
clear dcell;


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
mask_names = {'Viable Tumor Mask', 'Necrosis Mask', 'Vasculature Mask', 'Parenchyma Mask'};
outputVolumes = cell(4,3);
outputHeader = inputVolumes{1}{1};
outputHeader.RescaleSlope = 1;
outputHeader.RescaleIntercept = 0;
outputHeader.Type = 1; %1 for binary, else 0
outputMetaDataDict = {};
for i = 1:4
    outputVolumes{i}{1} = outputHeader;
    outputVolumes{i}{1}.Description = mask_names{i};
    outputVolumes{i}{2} = outputMetaDataDict;
    outputVolumes{i}{3} = masks{i};
end

status = true;
message = 'Successfully segmented input';

end

