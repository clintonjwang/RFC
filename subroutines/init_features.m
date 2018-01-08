function [ data ] = init_features( patients, working_dir )
%INIT_FEATURES Summary of this function goes here
%   Detailed explanation goes here

disp('Generating feature cell arrays (<1 min)...');
for i=1:length(patients)
    data{i} = load([working_dir,'/data_',num2str(i),'.mat']);
    data{i} = data{i}.data_i;
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
%     features{i} = generate_feature_array(data{i}); %returns only labels

    f = features{i};
    save([working_dir,'/init_features_',num2str(i),'.mat'], 'f');
end

end

