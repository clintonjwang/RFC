function compute_features(patients, working_dir)

if exist([working_dir,'/features_',num2str(1),'.mat'],'file') == 0
    tic
    for i=1:length(patients)
        data{i} = load([working_dir,'/data_',num2str(i),'.mat']);
        data{i} = data{i}.data_i;
        disp('Generating feature cell array...'); %fast
        features{i} = struct;
        features{i}.locations = find(data{i}.tight_liver_mask);
        features{i}.labels=zeros(length(features{i}.locations),1);

        for c=1:length(features{i}.locations)
            if(data{i}.necrosis_mask(features{i}.locations(c))==1)
                features{i}.labels(c)=3;
            elseif(data{i}.tumor_mask(features{i}.locations(c))==1)
                features{i}.labels(c)=1;
            elseif(data{i}.vessel_mask(features{i}.locations(c))==1)
                features{i}.labels(c)=2;
            end
        end
        features{i} = generate_feature_array(data{i}); %returns only labels

        f = features{i};
        save([working_dir,'/features_',num2str(i),'.mat'], 'f');
    end
    toc

    normalize_data(data, working_dir);  %fast
end

f = load([working_dir,'/features_',num2str(1),'.mat']);
f = f.f;
if ~isfield(f, 'auto_context_features_boost')
    for i=1:length(patients) %parfor
        compute_features_single(i, working_dir);
    end
end

return
end
