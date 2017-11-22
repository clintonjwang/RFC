function intensities_bin(patients, working_dir)
% train_size = 20;
% feat_dir = dir('./features');
% feat_dir = {feat_dir(3:end).name};
% feat_dir = ~cellfun(@(v) v(end-3:end) == '.mat', feat_dir, 'UniformOutput'=0);

% load([feat_dir,'/features_',num2str(i),'.mat']);

for i = 1:length(patients) %parfor
    f = load([working_dir,'/features_',num2str(i),'.mat']);
    f = f.f;
    
    intensities = reshape(f.intensities, [size(f.intensities,1)*size(f.intensities,2), 1]);
    
    fileID = fopen([working_dir,'/intensities_',num2str(i),'.bin'],'w');
    fwrite(fileID,intensities,'double');
    fclose(fileID);
    
    locations{i} = f.locations;
    f.num_intensity_features = size(f.intensities, 2); %301
    f.intensities = [];
    save([working_dir,'/features_',num2str(i),'.mat'],'f','-v7.3');
end
save([working_dir,'/locations.mat'],'locations','-v7.3');


return