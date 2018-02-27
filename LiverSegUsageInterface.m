function [ status, message, outputVolumes ] = LiverSegUsageInterface( pre, art, pv, t2, liver, model_dir, working_dir )
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

    %% Initialize
    status = false;

    data = struct;
    dims = pre{1}.NumberOfVoxels';
    data.pre = reshape(pre{3}, dims);
    data.art = reshape(art{3}, dims);
    data.pv = reshape(pv{3}, dims);
    data.t2 = reshape(t2{3}, dims);
    data.liver_seg = reshape(liver{3}, dims);
    
    filename_map = containers.Map;    
    filename_map('pre') = [working_dir, '/temp/pre_reg.nii'];
    filename_map('art') = [working_dir, '/temp/20s.nii'];
    filename_map('pv') = [working_dir, '/temp/70s_reg.nii'];
    filename_map('t2') = [working_dir, '/temp/t2_reg.nii'];
    filename_map('liver_seg') = [working_dir, '/temp/liver.nii'];
    for k = {'pre','art','pv','t2','liver_seg'}
        k = char(k);
        save_nii(data.(k), filename_map(k));
    end
    
    patients = {'0'};
    train_bool = false;

    %% Processing
    % get number of voxels and number of slices from the header
    acquire_data(patients, working_dir, working_dir, train_bool, filename_map);
    % Initialize labels and locations, and compute image features
    normalize_data(patients, working_dir);
    % Separate features based on label
    compute_features(patients, working_dir);
    % Save intensities in a separate bin file
    intensities_bin(patients, working_dir);
    % Train random forest model
    masks = tissue_classification(patients, model_dir, working_dir, train_bool, [], [], []);
    [~, ~, ~] = rmdir(working_dir, 's');

    mask_names = {'Viable Tumor Mask', 'Necrosis Mask', 'Vasculature Mask', 'Parenchyma Mask'};
    outputVolumes = cell(4,3);
    outputHeader = liver{1};
    for i = 1:4
        outputVolumes{i}{1} = outputHeader;
        outputVolumes{i}{1}.Description = mask_names{i};
        outputVolumes{i}{2} = outputMetaDataDict;
        outputVolumes{i}{3} = uint16(masks{i}(:));
    end

    status = true;
    message = 'Classification complete.';
end