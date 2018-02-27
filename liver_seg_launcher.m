function liver_seg_launcher
% Treilhard J. et al. (2017) Liver Tissue Classification in Patients with
% Hepatocellular Carcinoma by Fusing Structured and Rotationally Invariant
% Context Representation. In: Descoteaux M., Maier-Hein L., Franz A.,
% Jannin P., Collins D., Duchesne S. (eds) Medical Image Computing and
% Computer-Assisted Intervention - MICCAI 2017. MICCAI 2017. Lecture Notes
% in Computer Science, vol 10435. Springer, Cham
% DOI https://doi.org/10.1007/978-3-319-66179-7_10

    button = questdlg(['This program segments livers from T1/T2 MRIs '...
        'into viable tumor, necrosis, vasculature and parenchyma using '...
        'a cascading random forest model (Treilhard et al., MICCAI 2017). '...
        'Are you segmenting a batch of patients or a single patient?'], 'Liver Segmenter',...
        'Batch', 'Single patient',... %'Display mask',...
        'Use the trained model');
    
    switch button
        case 'Batch'
            user_main(button,false);
        case 'Single patient'
            user_main(button,false);
        case 'Display mask'
    %         disp('Not yet ready');
    %         mask_names = {'vasculature_mask', 'necrosis_mask', 'viable_tumor_mask'};
            mask_display_names = {'vasculature', 'necrosis', 'viable tumor'};
            mask_names = {'Vessels', 'necro', 'whole'};

            data_dir = uigetdir('', 'Select the folder containing the arterial image.');
            if data_dir == 0
                return
            end
            out_dir = uigetdir('', 'Select the folder containing the masks.');
            if out_dir == 0
                return
            end

            display_scrolling_mask('20s', data_dir, out_dir, mask_names, mask_display_names);

        case ''
            return
    end
end