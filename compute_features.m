function compute_features(patients, working_dir)

% f = load([working_dir,'/features_',num2str(1),'.mat']);
% f = f.f;

if exist([working_dir,'/features_',patients{1},'.mat'],'file') == 0
    disp('Computing features for each label...');
    tic
% if ~isfield(f, 'auto_context_features_boost')
    parfor i=1:length(patients)
        data = load([working_dir,'/norm_data_',patients{i},'.mat']);
        data = data.data_i;
        f = load([working_dir,'/init_features_',patients{i},'.mat']);
        f = f.f;
        f = compute_features_single(data, f);
%         save([working_dir,'/features_',patients{i},'.mat'],'f');
    end
    toc
end

return
end
