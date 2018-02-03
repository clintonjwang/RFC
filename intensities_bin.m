function intensities_bin(patients, working_dir)
% train_size = 20;
% feat_dir = dir('./features');
% feat_dir = {feat_dir(3:end).name};
% feat_dir = ~cellfun(@(v) v(end-3:end) == '.mat', feat_dir, 'UniformOutput'=0);

% load([feat_dir,'/features_',num2str(i),'.mat']);

disp('Saving intensities as bin files...');

locations = cell(1,length(patients));
parfor i = 1:length(patients)
    f = load_wrapper([working_dir,'/features_',patients{i},'.mat']);
    
    if exist([working_dir,'/intensities',patients{i},'.bin'],'file') == 0
        intensities = reshape(f.intensities, [size(f.intensities,1)*size(f.intensities,2), 1]);

        fileID = fopen([working_dir,'/intensities_',patients{i},'.bin'],'w');
        fwrite(fileID,intensities,'double');
        fclose(fileID);

        f.num_intensity_features = size(f.intensities, 2); %301
        f.intensities = [];
        save_wrapper(f, [working_dir,'/small_features_',patients{i},'.mat']);
    end
    
    locations{i} = f.locations;
end
save_wrapper(locations, [working_dir,'/locations.mat']);

return