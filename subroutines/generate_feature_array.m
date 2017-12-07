function features = generate_feature_array(data)
%generate_feature_array UNUSED

features = struct;

features.locations = find(data.tight_liver_mask);
features.labels=zeros(length(features.locations),1);

for c=1:length(features.locations)
    if(data.necrosis_mask(features.locations(c))==1)
        features.labels(c)=3;
    elseif(data.tumor_mask(features.locations(c))==1)
        features.labels(c)=1;
    elseif(data.vessel_mask(features.locations(c))==1)
        features.labels(c)=2;
    end
end

return
end