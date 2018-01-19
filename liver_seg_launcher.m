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
    'Are you retraining the model or using it?'], 'Liver Segmenter',...
    'Retrain the model (~30 min/patient)', 'Use the trained model (~30 min)', 'Display mask',...
    'Use the trained model (~30 min)');

switch button
    case 'Use the trained model (~30 min)'
        user_main(false);
    case 'Retrain the model (~30 min/patient)'
        train_main(false);
    case 'Display mask'
%         disp('Not yet ready');
%         mask_names = {'vasculature_mask', 'necrosis_mask', 'viable_tumor_mask'};
        mask_display_names = {'vasculature', 'necrosis', 'viable tumor'};
        mask_names = {'essels', 'necro', 'whole'};
        
        data_dir = 'E:/4-segmented lesions/ACcGBit/nii_files';
        out_dir = 'E:/4-segmented lesions/ACcGBit/new segs';
%         data_dir = uigetdir('', 'Select the folder containing the arterial image.');
%         if data_dir == 0
%             return
%         end
%         out_dir = uigetdir('', 'Select the folder containing the masks.');
%         if out_dir == 0
%             return
%         end
        
        display_scrolling_mask('20s', data_dir, out_dir, mask_names, mask_display_names);

    case ''
        return
end
end