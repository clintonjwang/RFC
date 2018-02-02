function f = compute_features_single(data, f, patient_id)
    num_intensity_maps = length(data.p_im);
    num_frangi_maps = length(data.frangi); 
    num_grad_maps = length(data.grad); 
    num_sf_maps = length(data.sf);
    num_mode_maps = length(data.mode);
    num_haralick = 3;

    f.sz = size(data.tight_liver_mask);

    % disp('Extracting voxel-wise intensities...'); %very fast
    f.intensities = zeros(length(f.locations),...
        num_intensity_maps);
    for j=1:num_intensity_maps
        f.intensities(:,j) = data.p_im{j}(f.locations); 
    end

    % disp('Extracting mode-map intensities...'); %very fast
    f.mode_intensities = zeros(length(f.locations),...
        num_mode_maps);
    for j=1:num_mode_maps
        f.mode_intensities(:,j) = data.mode{j}(f.locations);
    end

    % disp('Extracting Frangi intensities...');  %very fast
    f.frangi = zeros(length(f.locations),...
        num_frangi_maps);
    for j=1:num_frangi_maps
        f.frangi(:,j) = data.frangi{j}(f.locations); 
    end

    % disp('Extracting t2 intensities...'); %very fast
    f.t2 = zeros(length(f.locations),1);
    f.t2(:,1) = data.t2(f.locations);

    % disp('Extracting surface distances...'); %very fast
    f.intensities = [f.intensities,...
        compute_surface_distance(data,f.locations)];

    % disp('Extracting gradient intensities...');  %very fast
    f.gradient = zeros(length(f.locations),...
        num_grad_maps);
    for j=1:num_grad_maps
        f.gradient(:,j) = data.grad{j}(f.locations);
    end

    % disp('Extracting std filter intensities...'); %very fast
    f.sf = zeros(length(f.locations),...
        num_sf_maps);
    for j=1:num_sf_maps
        f.sf(:,j) = data.sf{j}(f.locations);
    end

    % disp('Extracting Haralick feature intensities...'); %very fast
    f.haralick = zeros(length(f.locations),num_intensity_maps*...
        num_haralick);

    ind=1;
    for j1=1:num_intensity_maps
        for j2=1:num_haralick
            f.haralick(:,ind) = data.glcm{j1,j2}(...
                f.locations);
            ind=ind+1;
        end
    end

    f.intensities = [f.intensities,...
        f.mode_intensities];
    f=rmfield(f,'mode_intensities');

    f.intensities = [f.intensities,...
        f.haralick];
    f=rmfield(f,'haralick');

    f.intensities = [f.intensities,...
        f.sf];
    f=rmfield(f,'sf');

    f.intensities = [f.intensities,...
        f.gradient];
    f=rmfield(f,'gradient');

    f.intensities = [f.intensities,...
        f.frangi];
    f=rmfield(f,'frangi');

    f.intensities = [f.intensities,...
        f.t2];
    f=rmfield(f,'t2');


    % disp('Appending contextual features'); %moderately fast
    f = append_context_features(f,data);

    f.auto_context_features=[]; 
    f.auto_context_features_boost=[]; 

    save_dummy(features, [save_dir,'/features_',patient_id,'.mat']);

    return
end

function [] =save_dummy(f,path)

save(path,'f','-v7.3');

return 
end

function features_i = append_context_features(features_i,data_i)

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