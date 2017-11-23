function compute_features(patients, working_dir)

% f = load([working_dir,'/features_',num2str(1),'.mat']);
% f = f.f;

if exist([working_dir,'/features_',num2str(1),'.mat'],'file') == 0
    disp('Computing features for each label...');
    tic
% if ~isfield(f, 'auto_context_features_boost')
    for i=1:length(patients) %parfor
        compute_features_single(i, working_dir);
    end
    toc
end

return
end
