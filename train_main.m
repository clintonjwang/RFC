function train_main(skipgui)
%TRAIN_MAIN entry point for training the random forest tissue classifier.

    if nargin < 1
        skipgui = true;
    end

    addpath(genpath('../subroutines'));
    
    train_bool = true;
    model_dir = fullfile(pwd(),'models');
    working_dir = 'D:/working_train';

    filename_map = containers.Map;
    filename_map('pre') = '**/pre_reg.nii*';
    filename_map('art') = '**/20s.nii*';
    filename_map('pv') = '**/70s_reg.nii*';
    filename_map('t2') = '**/t2_bfc_reg.nii*';
    filename_map('t1-bfc') = 'nii_files/bias_field_isotropic.nii';
    filename_map('liver_seg') = '**/*liver.ids';
    filename_map('tumor_seg') = '**/*tumor*.ids';
    filename_map('vasc_seg') = '**/*vessel*.ids';
    filename_map('necro_seg') = '**/*nec*.ids';

    %% Take user input
    % Instructions
    if ~skipgui
        uiwait(msgbox(['All training data should be stored in a single folder '...
            'where each subfolder contains the data for a single patient. '...
            'There should be no extraneous subfolders. '...
            'Be sure that each patient subfolder has axial precontrast, arterial, '...
            'and venous T1-w images, axial T2-w images, a binary mask '...
            'representing the whole liver, as well as 4 binary masks '...
            'representing whole liver, tumor, necrosis, and vasculature '...
            'segmentations, each in ICS version 1 format.'...
            'Once started, it may take over an hour per patient.'...
            'WARNING: At its peak, it uses around 7GB of memory.'], 'Random Forest training', 'modal'));
    end
    
    % Set training directory
    if skipgui
        data_dir = 'E:/4-segmented HCCs';
    else
        data_dir = uigetdir('', 'Select the folder containing all the training data.');
        if data_dir == 0
            return
        end
    end

    if false%~skipgui
        prompt = {'Save features at the end of the run?',...
                'Do images have separate T1-w bias field corrections saved as a nifti?'};
        dlg_title = 'Run options';
        num_lines = 1;
        defaultans = {'no','yes'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

        if ~iscell(answer)
            return
        end
    else
        answer = {'yes','yes'};
    end

    save_features = strcmp(answer{1},'yes') == 1;
    use_bias_field = strcmp(answer{2},'yes') == 1;

    % Specify training options
    if ~skipgui
        prompt = {'Enter number of decision trees in each random forest',...
                'Enter number of classification rounds',...
                'Enter structured context patch size',...
                'Enter spherical context patch size',...
                'Enter number of histogram bins for spherical histogram context features',...
                'Enter number of auto-context features to sample per dimension',...
                'Enter mininum number of leaves in decision trees'};
        dlg_title = 'Cascading random forest parameters';
        num_lines = 1;
        defaultans = {'800','3','8','5','5','6','50'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

        if ~iscell(answer)
            return
        end
    else
        answer = {'800','3','8','5','5','6','50'};
    end

    config.ntrees = str2double(answer{1});
    config.RAC = str2double(answer{2});
    config.sl = str2double(answer{3});
    config.sl_spherical = str2double(answer{4});
    config.num_bins = str2double(answer{5});
    config.sa = str2double(answer{6});
    config.min_leaf_size = str2double(answer{7});

    [~, ~, ~] = mkdir(model_dir);
    [~, ~, ~] = mkdir(working_dir);

    patients = dir(data_dir);
    filenames = {patients.name};
    patients = filenames([patients.isdir]);
    patients = patients(3:end);

    %% Run algorithm
    tic
    % Collect images and whole liver masks
    acquire_data(patients, data_dir, working_dir, train_bool, use_bias_field, filename_map);
    toc

    % Initialize labels and locations, and compute image features
    normalize_data(patients, working_dir);
    toc

    % Separate features based on label
    compute_features(patients, working_dir);
    toc

    % Generate training data
    if exist([working_dir,'/voxel_data.mat'],'file') == 0
        generate_training_data(patients, working_dir);
        toc
    end

    % Save intensities in a separate bin file
    intensities_bin(patients, working_dir);
    toc

    % Train random forest model
    tissue_classification(patients, model_dir, working_dir, train_bool, data_dir, [], config);
    toc

    if ~save_features
        [~, ~, ~] = rmdir(working_dir, 's');
    end

    if ~skipgui
        uiwait(msgbox(['Training complete. The tree weights have been saved to ',...
            model_dir, '.'], 'Random Forest training complete', 'modal'));
    end

end