feat_dir = dir('/gpfs/home/fas/duncan/jy498/Feature');
feat_dir = {feat_dir(3:57).name};

for i = 1:length(feat_dir)
    i
    load(['/gpfs/home/fas/duncan/jy498/Feature/features_',num2str(i),'.mat']);
    intensities = reshape(f.intensities,[size(f.intensities,1)*size(f.intensities,2),1]);
    fileID = fopen(['/gpfs/scratch60/fas/duncan/jy498/Feature/intensities_',num2str(i),'.bin'],'w');
    fwrite(fileID,intensities,'double');
    fclose(fileID);
    locations{i} = f.locations;
    f.num_intensity_features = size(f.intensities,2);
    f.intensities = [];
    save(['/gpfs/home/fas/duncan/jy498/Research/Junlin/Feature/features_',num2str(i),'.mat'],'f','-v7.3');

end   

save('/gpfs/home/fas/duncan/jy498/Research/Junlin/Parameters/locations.mat','locations','-v7.3');