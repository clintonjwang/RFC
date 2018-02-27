function user_main(mode, skipgui)
%USER_MAIN entry point for using the trained random forest classifier
    
    if ~exist('button','var') || isempty(mode)
        mode = 'Single patient';
    end
    if ~exist('skipgui','var') || isempty(skipgui)
        skipgui = true;
    end

    addpath(genpath('../subroutines'));

    train_bool = false;
    model_dir = fullfile(pwd(),'models');
    mask_dir = fullfile(pwd(),'masks');
    working_dir = fullfile(pwd(),'working_test');

    %% Take user input
%     if ~skipgui
%         uiwait(msgbox(['Using this model requires that it has been previously '...
%             'trained. This info should be saved as "tree_x.mat". '...
%             'The model also requires registered contrast enhanced T1/T2 '...
%             'MRIs (in nifti format) along with the whole liver mask (in .ics/'...
%             '.ids format). It outputs binary masks of the 4 classes (also in '...
%             '.ics/.ids format). Only ICS version 1 is supported. '...
%             'The program asks you to select a patient folder to look for the '...
%             'MRIs/masks in. If it cannot find a file automatically, it will '...
%             'prompt you for it.'], 'Random Forest Cascade utility', 'modal'));
%     end
    if strcmp(mode, 'Single patient')
        config = load_config(try_find_file(pwd(), '**/config.txt', ...
            'Select the config file.', {'.txt'}));
        
        filename_map = containers.Map;
        for k = {'pre','art','pv','t2','liver_seg'}
            k = char(k);
            filename_map(k) = config(k);
        end
        data_dir = working_dir;
        model_dir = config('models');
        mask_dir = config('mask_output');
        patients = {'0'};
        
    else
        if skipgui
            data_dir = 'E:/4-sampler';
        else
            data_dir = uigetdir('', ['Select the folder containing all patients to segment. '...
                    '(Choose a parent folder containing each patient as a separate subfolder.)']);
            if data_dir == 0
                return
            end
        end

        if ~exist(model_dir,'dir')
            model_dir = uigetdir('', 'Select the folder containing all the models.');
            if model_dir == 0
                return
            end
        end

        % Specify file naming convention
        defaultans = {'pre_reg.nii*','20s.nii*','70s_reg.nii*','t2_bfc_reg.nii*','*liver.ids'};
        if ~skipgui
            prompt = {'Enter file pattern for pre-contrast nifti',...
                    'Enter file pattern for arterial phase nifti',...
                    'Enter file pattern for venous phase nifti',...
                    'Enter file pattern for T2-w nifti',...
                    'Enter file pattern for liver mask'};
            dlg_title = 'File naming convention';
            num_lines = 1;
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

            if ~iscell(answer)
                return
            end
        else
            answer = defaultans;
        end
        filename_map = containers.Map;
        filename_map('pre') = ['**/', answer{1}];
        filename_map('art') = ['**/', answer{2}];
        filename_map('pv') = ['**/', answer{3}];
        filename_map('t2') = ['**/', answer{4}];
        filename_map('liver_seg') = ['**/', answer{5}];

        if skipgui
            mask_dir = fullfile(pwd(),'masks');
        else
            mask_dir = uigetdir('', 'Select a folder to output the segmentations to.');
            if mask_dir == 0
                return
            end
        end

        patients = dir(data_dir);
        filenames = {patients.name};
        patients = filenames([patients.isdir]);
        patients = patients(3:end);
    end

    defaultans = {'no'};
    if ~skipgui
        prompt = {'Save features at the end of the run?'};%,...
                %'Is there a separate T1-w bias field correction?'};
        dlg_title = 'Run options';
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

        if ~iscell(answer)
            return
        end
    else
        answer = defaultans;
    end

    save_features = strcmp(answer{1},'yes') == 1;

    [~, ~, ~] = mkdir(mask_dir);
    [~, ~, ~] = mkdir(working_dir);

    %% Run algorithm
    tic
    % Collect images and whole liver masks
    acquire_data(patients, data_dir, working_dir, train_bool, filename_map);
    toc

    % Initialize labels and locations, and compute image features
    normalize_data(patients, working_dir);
    toc

    % Separate features based on label
    compute_features(patients, working_dir);
    toc

    % Save intensities in a separate bin file
    intensities_bin(patients, working_dir);
    toc

    % Train random forest model
    tissue_classification(patients, model_dir, working_dir, train_bool, data_dir, mask_dir, []);
    toc

    if ~save_features
        [~, ~, ~] = rmdir(working_dir, 's');
    end

    if ~skipgui
        uiwait(msgbox(['Classification complete. The masks have been saved to ',...
            mask_dir, '.'], 'Random Forest complete', 'modal'));
    end

    % Display result
%     mask_names = {'vasculature_mask', 'necrosis_mask', 'viable_tumor_mask'};
%     mask_display_names = {'vasculature', 'necrosis', 'viable tumor'};
%     display_scrolling_mask('20s', data_dir, out_dir, mask_names, mask_display_names);
end