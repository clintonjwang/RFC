function normalize_data(data, working_dir)
%normalize_data(data,features) computes features for classification 
% data should contain all patients

disp('Normalizing data...');

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
for i=1:length(data) %parfor
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
    
    disp('Computing gradient and std filter maps...'); %fast (<1s)
    for j=1:length(data{i}.p_im)
        data{i}.grad{j} = compute_gradient_image(data{i}.p_im{j});
        data{i}.sf{j} = stdfilt(data{i}.p_im{j});
    end
    
    disp('Computing Frangi features...'); %fast (<1 min)
    tic
    data{i}.frangi{1} = FrangiFilter3D(data{i}.p_im{3},Options1);
    data{i}.frangi{2} = FrangiFilter3D(data{i}.p_im{2},Options2);
    data{i}.frangi{3} = FrangiFilter3D(data{i}.p_im{1},Options1);
    
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

disp('Computing texture features...');  % slow
for i=1:length(data) %parfor
    tic
    f = load([working_dir,'/features_',num2str(i),'.mat']);
    f = f.f;
    data{i}.glcm = compute_glcm(data{i}.p_im,f.locations,range); 
    toc
    
    data_i = data{i};
    save([working_dir,'/data_',num2str(i),'.mat'], 'data_i');
end

return
end


function gim = compute_gradient_image(im) 
%gradient image

[FX,FY,FZ]=gradient(im); 
gim = ((FX.^2)+(FY.^2)+(FZ.^2)).^0.5;

return 
end

function glcm_maps = compute_glcm(maps,locations,range)
%haralick texture feature image

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


%compute image where intensity corresponds to distance from image mode 
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