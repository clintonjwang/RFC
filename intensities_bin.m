feat_dir = dir('./features');
feat_dir = {feat_dir(3:end).name};

for i = 1:length(feat_dir)
    i
    load(['./features/features_',num2str(i),'.mat']);
    intensities = reshape(f.intensities,[size(f.intensities,1)*size(f.intensities,2),1]);
    fileID = fopen(['./features/intensities_',num2str(i),'.bin'],'w');
    fwrite(fileID,intensities,'double');
    fclose(fileID);
    locations{i} = f.locations;
    f.num_intensity_features = size(f.intensities,2);
    f.intensities = [];
    save(['./features/features_',num2str(i),'.mat'],'f','-v7.3');

end   

save('./params/locations.mat','locations','-v7.3');