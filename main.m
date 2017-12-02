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
    'Retrain the model', 'Use the trained model', 0);

switch button
    case 'Use the trained model'
        user_main(false)
    case 'Retrain the model'
        train_main(false)
        disp('k')
    case 0
        disp('eagargarg')
        return
end