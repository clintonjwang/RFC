function intensities_bin(f)
% train_size = 20;
feat_dir = 'features';
% feat_dir = dir('./features');
% feat_dir = {feat_dir(3:end).name};
% feat_dir = ~cellfun(@(v) v(end-3:end) == '.mat', feat_dir, 'UniformOutput'=0);

% load([feat_dir,'/features_',num2str(i),'.mat']);

intensities = reshape(f.intensities, [size(f.intensities,1)*size(f.intensities,2), 1]);
f.num_intensity_features = size(f.intensities, 2);
f.intensities = [];

return