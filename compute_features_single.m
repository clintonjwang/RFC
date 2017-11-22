function compute_features_single(data_path, pr, features)

data = load(data_path);
data = data.data_i;
num_intensity_maps = length(data.p_im);
num_frangi_maps = length(data.frangi); 
num_grad_maps = length(data.grad); 
num_sf_maps = length(data.sf);
num_mode_maps = length(data.mode);
num_haralick = 3;

features.sz = size(data.tight_liver_mask);

tic
disp('Extracting voxel-wise intensities...'); %very fast
features.intensities = zeros(length(features.locations),...
    num_intensity_maps);
for j=1:num_intensity_maps
    features.intensities(:,j) = data.p_im{j}(features.locations); 
end

disp('Extracting mode-map intensities...'); %very fast
features.mode_intensities = zeros(length(features.locations),...
    num_mode_maps);
for j=1:num_mode_maps
    features.mode_intensities(:,j) = data.mode{j}(features.locations);
end

disp('Extracting Frangi intensities...');  %very fast
features.frangi = zeros(length(features.locations),...
    num_frangi_maps);
for j=1:num_frangi_maps
    features.frangi(:,j) = data.frangi{j}(features.locations); 
end

disp('Extracting t2 intensities...'); %very fast
features.t2 = zeros(length(features.locations),1);
features.t2(:,1) = data.t2(features.locations);

disp('Extracting surface distances...'); %very fast
features.intensities = [features.intensities,...
    compute_surface_distance(data,features.locations)];

disp('Extracting gradient intensities...');  %very fast
features.gradient = zeros(length(features.locations),...
    num_grad_maps);
for j=1:num_grad_maps
    features.gradient(:,j) = data.grad{j}(features.locations);
end

disp('Extracting std filter intensities...'); %very fast
features.sf = zeros(length(features.locations),...
    num_sf_maps);
for j=1:num_sf_maps
    features.sf(:,j) = data.sf{j}(features.locations);
end

disp('Extracting Haralick feature intensities...'); %very fast
features.haralick = zeros(length(features.locations),num_intensity_maps*...
    num_haralick);

ind=1;
for j1=1:num_intensity_maps
    for j2=1:num_haralick
        features.haralick(:,ind) = data.glcm{j1,j2}(...
            features.locations);
        ind=ind+1;
    end
end

features.intensities = [features.intensities,...
    features.mode_intensities];
features=rmfield(features,'mode_intensities');

features.intensities = [features.intensities,...
    features.haralick];
features=rmfield(features,'haralick');

features.intensities = [features.intensities,...
    features.sf];
features=rmfield(features,'sf');

features.intensities = [features.intensities,...
    features.gradient];
features=rmfield(features,'gradient');

features.intensities = [features.intensities,...
    features.frangi];
features=rmfield(features,'frangi');

features.intensities = [features.intensities,...
    features.t2];
features=rmfield(features,'t2');
toc


disp('Appending contextual features');
tic
features = append_context_features(features,data,pr);
toc

features.auto_context_features=[]; 
features.auto_context_features_boost=[]; 

save_dummy(features,1);

features=features.labels;

return
end

function [] =save_dummy(f,i)

save_dest='features';
save([save_dest,'/features_',num2str(i),'.mat'],'f','-v7.3');

return 
end

function features_i = append_context_features(features_i,data_i,pr)

maps{1}=data_i.mode{1}; 
maps{2}=data_i.mode{2};
maps{3}=data_i.mode{3};
maps{4}=data_i.mode{4};

sl=4;

features_i.intensities = [features_i.intensities,...
    compute_patches(maps,features_i,sl,data_i)];

return 
end

function temp = compute_patches(maps,features_i,sl,data_i)
[ii,jj,kk] = ind2sub(features_i.sz,features_i.locations);
num_maps = length(maps);
N=4;

for c=1:length(maps)
    patches{c}=zeros(length(ii),N^3); 
end

for d=1:length(ii)
    
    x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,features_i.sz(1)),N));
    y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,features_i.sz(2)),N));
    z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,features_i.sz(3)),N));
    sample_locs = combvec(x,y,z);
    
    for n=1:(N^3)
        if(data_i.tight_liver_mask(sample_locs(1,n),sample_locs(2,n),...
                sample_locs(3,n))==0)
            sample_locs(:,n)=[ii(d);jj(d);kk(d)];
        end
    end
    
    inds = sub2ind(features_i.sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));

    for c=1:length(maps)
        patches{c}(d,:) = maps{c}(inds); 
    end
    
end

temp=[];
for c=1:num_maps
   temp=[temp,patches{c}];  
end

return
end


function temp = compute_surface_distance(data_i,locations)

D=bwdist(1-data_i.tight_liver_mask);

temp=zeros(length(locations),1);
for j=1:length(locations)
    temp(j)=D(locations(j));
end

temp = reshape(temp,[numel(temp),1]); 

return
end

%{
function temp = compute_patches_spherical_distributions(...
    maps,features_i,sl,data_i,pr)

M = generate_spherical_masks(sl);
num_spheres=length(M);
[ii,jj,kk] = ind2sub(features_i.sz,features_i.locations);
num_maps = length(maps);
N=(2*sl)+1;
num_bins=25; 

for c=1:length(maps)
    edges{c} = [-Inf,linspace(pr{c}(1),pr{c}(2),num_bins-1),+Inf];
end

for c=1:length(maps)
    patches{c}=zeros(length(ii),num_bins*num_spheres); 
end

for d=1:length(ii)
    
    x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,features_i.sz(1)),N));
    y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,features_i.sz(2)),N));
    z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,features_i.sz(3)),N));
    sample_locs = combvec(x,y,z);
    
    for n=1:(N^3)
        if(data_i.tight_liver_mask(sample_locs(1,n),sample_locs(2,n),...
                sample_locs(3,n))==0)
            sample_locs(:,n)=[ii(d);jj(d);kk(d)];
        end
    end
    
    inds = sub2ind(features_i.sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));
    
    for c=1:length(maps)
        for s=1:num_spheres
            patches{c}(d,((s-1)*num_bins)+1:s*num_bins) = ...
                histcounts(maps{c}(inds(M{s})),edges{c});
        end
    end
    
end

temp=[];
for c=1:num_maps
   temp=[temp,patches{c}];  
end

return
end
%}