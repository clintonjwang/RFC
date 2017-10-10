%compute features for classification 
function pr = normalize_data(data,features)

%parameters for Frangi vessel filter computation 
Options1.FrangiScaleRange=[1,10];
Options1.FrangiScaleRatio=1;
Options1.verbose=0;
Options1.BlackWhite=0;

Options2.FrangiScaleRange=[1,10];
Options2.FrangiScaleRatio=1;
Options2.verbose=0;
Options2.BlackWhite=1;

%loop through patients 
parfor i=1:length(data)
%for i=1:length(data)
    tic
    
    %bias field correction 
    data{i}.pre = data{i}.pre./data{i}.bf;
    data{i}.art = data{i}.art./data{i}.bf;
    data{i}.pv = data{i}.pv./data{i}.bf; 
    
    %normalize data by making it 0 mean, unit variance 
    m = mean(data{i}.pre(data{i}.tight_liver_mask==1));
    s = sqrt(var(data{i}.pre(data{i}.tight_liver_mask==1)));
    
    data{i}.pre = data{i}.pre - m;
    data{i}.art = data{i}.art - m;
    data{i}.pv = data{i}.pv - m;
    
    data{i}.pre = data{i}.pre / s;
    data{i}.art = data{i}.art / s;
    data{i}.pv = data{i}.pv / s;
    
    data{i}.pre = data{i}.pre + 10;
    data{i}.art = data{i}.art + 10;
    data{i}.pv = data{i}.pv + 10;
    
    data{i}.p_im{1} = (data{i}.art - data{i}.pre)./data{i}.pre;
    data{i}.p_im{2} = (data{i}.art - data{i}.pv) ./ data{i}.art;
    data{i}.p_im{3} = (data{i}.pv - data{i}.pre) ./ data{i}.pre;
    data{i}.p_im{4} = data{i}.pre;
    data{i}.p_im{5} = data{i}.art;
    data{i}.p_im{6} = data{i}.pv;
    
    data{i}.mode{1} = compute_mode_map(data{i}.p_im{1},data{i}.tight_liver_mask);
    data{i}.mode{2} = compute_mode_map(data{i}.p_im{2},data{i}.tight_liver_mask);
    data{i}.mode{3} = compute_mode_map(data{i}.p_im{3},data{i}.tight_liver_mask);
    data{i}.mode{4} = compute_mode_map(data{i}.t2,data{i}.tight_liver_mask);
    
    %normalize T2 image to be 0 mean, unit variance 
    m = mean(data{i}.t2(data{i}.tight_liver_mask==1));
    s = sqrt(var(data{i}.t2(data{i}.tight_liver_mask==1)));
    data{i}.t2 = data{i}.t2 - m;
    data{i}.t2 = data{i}.t2 / s; 
    
    disp('Computing gradient and std filter maps...'); 
    tic
    for j=1:length(data{i}.p_im)
        data{i}.grad{j} = compute_gradient_image(data{i}.p_im{j});
        data{i}.sf{j} = stdfilt(data{i}.p_im{j});
    end
    toc
    
    disp('Computing Frangi features...'); 
    tic
    data{i}.frangi{1} = FrangiFilter3D(data{i}.p_im{3},Options1);
    data{i}.frangi{2} = FrangiFilter3D(data{i}.p_im{2},Options2);
    data{i}.frangi{3} = FrangiFilter3D(data{i}.p_im{1},Options1);
    toc
    
    if(isfield(data{i},'pre'))
        data{i} = rmfield(data{i},'pre');
    end
    
    if(isfield(data{i},'art'))
        data{i} = rmfield(data{i},'art');
    end
    
    if(isfield(data{i},'pv'))
        data{i} = rmfield(data{i},'pv');
    end
    
     if(isfield(data{i},'bf'))
        data{i} = rmfield(data{i},'bf');
    end
    toc
    
end

for j=1:length(data{1}.p_im)
    temp{j}=[];
end

for i=1:length(data)
    for j=1:length(data{1}.p_im)
        temp{j} = [temp{j};...
            data{i}.p_im{j}(data{i}.tight_liver_mask==1)];
    end
end

for j=1:length(data{1}.p_im)
   range{j}=prctile(temp{j},[1,99]); 
end
clear temp

disp('Computing texture features...'); 
parfor i=1:length(data)
    i
    tic
    data{i}.glcm = compute_glcm(data{i}.p_im,features{i}.locations,range); 
    toc
end

%extract 1st and 9th percentile of data
for i=1:3
    maps{i}=[]; 
end

for i=1:length(data) 
    maps{1}=[maps{1};data{i}.p_im{1}(:)];
    maps{2}=[maps{2};data{i}.p_im{2}(:)];
    maps{3}=[maps{3};data{i}.p_im{3}(:)];
end

for i=1:3
    pr{i} = prctile(maps{i},[1,99]);
    minmax{i}=[min(maps{i}),max(maps{i})];
end

clear maps

%place to save data files after processing 
save_dir='/gpfs/home/fas/duncan/jy498/Research/Junlin/Data';

for i=1:length(data)
    data_i = data{i};
    save([save_dir,'/data_',num2str(i),'.mat'],'data_i','-v7.3'); 
    data{i}=[];
end

return
end

%gradient image 
function gim = compute_gradient_image(im) 

[FX,FY,FZ]=gradient(im); 
gim = ((FX.^2)+(FY.^2)+(FZ.^2)).^0.5;

return 
end

%haralick texture feature image 
function glcm_maps = compute_glcm(maps,locations,range)

num_haralick=3; 
sl=2;
sn=(2*sl)+1;
sz=size(maps{1}); 
[ii,jj,kk] = ind2sub(sz,locations);
num_maps = length(maps);

for c=1:num_maps
    for h=1:num_haralick
        glcm_maps{c,h} = zeros(size(maps{1}));
    end
end

for d=1:length(ii)
    
    x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,sz(1)),sn)); 
    y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,sz(2)),sn));
    
    for c=1:length(maps)
        stats = graycoprops(graycomatrix(maps{c}(x,y,kk(d)),'GrayLimits',...
            range{c},'NumLevels',25));  
        glcm_maps{c,1}(ii(d),jj(d),kk(d)) = stats.Contrast;
        glcm_maps{c,2}(ii(d),jj(d),kk(d)) = stats.Energy;
        glcm_maps{c,3}(ii(d),jj(d),kk(d)) = stats.Homogeneity;
    end
    
end

return
end

%smooth image
function data_i = smooth_data(data_i)

num_iter=3;
delta_t=0.0682;
kappa=30;
option=1;

m = mean(data_i.pv(data_i.tight_liver_mask==1));

phases={'pre','art','pv','t2'};

voxel_spacing = data_i.resolution;
for j=1:length(phases)
    m1 = mean(data_i.(phases{j})(data_i.tight_liver_mask==1));
    data_i.(phases{j}) = data_i.(phases{j}) * (m/m1);
    data_i.(phases{j})=(m1/m)*anisodiff3D(data_i.(phases{j}) * (m/m1),...
        num_iter,delta_t,kappa,option,voxel_spacing);
end

return
end

%compute image where intensity corresponds from distance from image mode 


function mode_map = compute_mode_map(map,mask) 

%set up histogram, find the mode 
num_bins=200;
intensities = map(mask==1); 
prc = prctile(intensities,[1,99]); 
edges = [-Inf,linspace(prc(1),prc(2),num_bins-1),+Inf];
edge_dist = edges(3)-edges(2); 
[N,~]=histcounts(intensities,edges); 
[~,mode_bin]=max(N); 

%set up mode map 
mode_map=zeros(size(map)); 

%fill in mode map 
locations=find(mask); 
for p=locations
    mode_map(p) = compute_mode_distance(map(p),edges,mode_bin,...
        num_bins,edge_dist); 
end

return
end
 
function md = compute_mode_distance(value,edges,mode_bin,num_bins,...
    edge_dist)

if(value<=edges(2))
    md = 1-mode_bin;
elseif(value>=edges(end-1))
    md = num_bins - mode_bin; 
else
    md = (floor((value - edges(2))/edge_dist) + 2) - mode_bin; 
end

return
end