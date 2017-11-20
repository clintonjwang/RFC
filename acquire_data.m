function data = acquire_data(patients, path, train_bool)
%acquire_data(patients, path, train_bool)
% patients should be cell array of subfolder names each containing a patient
% path is the path to the patient folders
% set train_bool to true if training

disp('Acquiring patient data...');

%parfor i=1:length(patients)
for i=1:length(patients)
    data{i} = acquire_data_single_pat([path,'/',patients{i}], train_bool);
    disp([patients{i}, ' acquired!']);
end

return
