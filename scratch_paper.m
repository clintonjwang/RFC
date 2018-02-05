tic
parfor td=1:num_patients
    %load corresponding patient
    f = load_wrapper([working_dir,'/small_features_',patients{td},'.mat']);

    %read in appearance feature binary file for patient 
    fileID = fopen([working_dir, '/intensities_',patients{td},'.bin'],'r');
    intensities = fread(fileID,[numel(f.labels),...
        f.num_intensity_features],'double');
    fclose(fileID);
end
toc
tic
parfor td=1:num_patients
    f = load_wrapper([working_dir,'/features_',patients{td},'.mat']);
end
toc