function acquire_data(patients, train_dir, working_dir, train_bool)
%ACQUIRE_DATA(patients, path, train_bool)
% patients should be cell array of subfolder names each containing a patient
% path is the path to the patient folders
% set train_bool to true if training

disp('Acquiring patient data...');
%parfor i=1:length(patients)
for i=1:length(patients)
    if exist([working_dir,'/data_',num2str(i),'.mat'],'file') == 0
        data_i = acquire_data_single_pat([train_dir,'/',patients{i}], train_bool);
        save([working_dir,'/data_',num2str(i),'.mat'],'data_i')
        disp([patients{i}, ' acquired!']);
    end
end

return
