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
    'Retrain the model', 'Use the trained model', 'Display mask',...
    'Use the trained model');

switch button
    case 'Use the trained model'
        user_main(false);
    case 'Retrain the model'
        train_main(false);
    case 'Display mask'
        fig = uifigure;
        cbx = uicheckbox(fig, 'Text','Show Value',...
                          'Value', 1,...
                          'Position',[150 50 102 15]);
    case ''
        return
end
end