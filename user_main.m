function masks = user_main(skipgui)
%USER_MAIN entry point for using the trained random forest classifier

    if nargin < 1
        skipgui = true;
    end

    addpath(genpath('subroutines'));

    filename_map = containers.Map;
    filename_map('pre') = '**/pre_reg.nii*';
    filename_map('art') = '**/20s.nii*';
    filename_map('pv') = '**/70s_reg.nii*';
    filename_map('t2') = '**/t2_bfc_reg.nii*';
    filename_map('t1-bfc') = '**/bias_field_isotropic.nii';
    filename_map('liver_seg') = '**/*liver.ids';
    filename_map('tumor_seg') = '**/*whole tumor.ids';
    filename_map('vasc_seg') = '**/*vessel*.ids';
    filename_map('necro_seg') = '**/*nec*.ids';

    %% Take user input
    if ~skipgui
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

    if skipgui
        data_dir = 'E:/4-segmented HCCs/~a90ede8a';
    else
        data_dir = uigetdir('', 'Select the folder containing the patient to segment.');
        if data_dir == 0
            return
        end
    end

    train_bool = false;
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

    if skipgui
        out_dir = fullfile(data_dir,'output_masks');
        [~, ~, ~] = mkdir(out_dir);
    else
        out_dir = uigetdir('', 'Select a folder to output the binary masks to.');
        if out_dir == 0
            return
        end
    end

    if ~skipgui
        prompt = {'Save features at the end of the run?',...
                'Is there a separate T1-w bias field correction?'};
        dlg_title = 'Run options';
        num_lines = 1;
        defaultans = {'no','no'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

        if ~iscell(answer)
            return
        end
    else
        answer = {'yes','no'};
    end

    save_features = strcmp(answer{1},'yes') == 1;
    use_bias_field = strcmp(answer{2},'no') == 1;


    %% Run algorithm
    % Collect images and whole liver masks
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
    if exist([working_dir,'/train.mat'],'file') == 0
        generate_training_data(patients, working_dir);
        toc
    end

    % Save intensities in a separate bin file
    intensities_bin(patients, working_dir);
    toc
    
    % Train random forest model
    masks = tissue_classification({''}, model_dir, working_dir, working_dir, train_bool, data_dir, out_dir);

    if ~save_features
        [~, ~, ~] = rmdir(working_dir, 's');
    end

    % Display result
    mask_names = {'vasculature_mask', 'necrosis_mask', 'viable_tumor_mask'};
    mask_display_names = {'vasculature', 'necrosis', 'viable tumor'};
    display_scrolling_mask('20s', data_dir, out_dir, mask_names, mask_display_names);
end