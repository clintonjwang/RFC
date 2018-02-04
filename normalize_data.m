function normalize_data(patients, working_dir)
    %normalize_data(data,features) computes features for classification 
    % data should contain all patients
    % data gains p_im, which captures pre, art, pv intensities and their diffs
    % data gains mode, grad, sf, frangi, and glcm (Haralick)
        
    disp('Normalizing data...');
    
    skip = true;
    for i=1:length(patients)
        if ~exist([working_dir,'/norm_data_',patients{i},'.mat'],'file')
            skip = false;
            break
        end
    end
    if skip
        return
    end

    %parameters for Frangi vessel filter computation 
    Options1.FrangiScaleRange=[1,10];
    Options1.FrangiScaleRatio=1;
    Options1.verbose=0;
    Options1.BlackWhite=0;

    Options2.FrangiScaleRange=[1,10];
    Options2.FrangiScaleRatio=1;
    Options2.verbose=0;
    Options2.BlackWhite=1;

    % create feature cells for each class type
   init_features(patients, working_dir);
    
    %loop through patients
    parfor i=1:length(patients)
        if ~exist([working_dir,'/norm_data_',patients{i},'.mat'],'file')
            data = load_wrapper([working_dir,'/data_',patients{i},'.mat']);
            
            if ~isfield(data, 'pre')
                continue
            end
            %bias field correction
            if isfield(data, 't1_bfc')
                data.pre = data.pre./data.t1_bfc;
                data.art = data.art./data.t1_bfc;
                data.pv = data.pv./data.t1_bfc;
            end

            %normalize data by making it 0 mean, unit variance 
            m = mean(data.pre(data.tight_liver_mask==1));
            s = sqrt(var(data.pre(data.tight_liver_mask==1)));

            for field = {'pre', 'art', 'pv'}
                data.(char(field)) = (data.(char(field)) - m) / s + 10;
            end

            data.p_im{1} = (data.art - data.pre)./data.pre;
            data.p_im{2} = (data.art - data.pv) ./ data.art;
            data.p_im{3} = (data.pv - data.pre) ./ data.pre;
            data.p_im{4} = data.pre;
            data.p_im{5} = data.art;
            data.p_im{6} = data.pv;

            data.mode{1} = compute_mode_map(data.p_im{1},data.tight_liver_mask);
            data.mode{2} = compute_mode_map(data.p_im{2},data.tight_liver_mask);
            data.mode{3} = compute_mode_map(data.p_im{3},data.tight_liver_mask);
            data.mode{4} = compute_mode_map(data.t2,data.tight_liver_mask);

            m = mean(data.t2(data.tight_liver_mask==1));
            s = sqrt(var(data.t2(data.tight_liver_mask==1)));
            data.t2 = (data.t2 - m) / s;

            % Computing gradient and std filter maps
            for j=1:length(data.p_im)
                data.grad{j} = compute_gradient_image(data.p_im{j});
                data.sf{j} = stdfilt(data.p_im{j});
            end

            % Computing Frangi features
            data.frangi{1} = FrangiFilter3D(data.p_im{3},Options1);
            data.frangi{2} = FrangiFilter3D(data.p_im{2},Options2);
            data.frangi{3} = FrangiFilter3D(data.p_im{1},Options1);
            
            for k = {'pre','art','pv','t1_bfc'}
                k = char(k);
                if(isfield(data,k))
                    data = rmfield(data,k);
                end
            end
            save_wrapper(data, [working_dir,'/data_',patients{i},'.mat']);
        end
    end

    
    data = load_wrapper([working_dir,'/data_',patients{1},'.mat']);
    for j=1:length(data.p_im)
        temp{j}=[];
    end

    for i=1:length(patients)
        data = load_wrapper([working_dir,'/data_',patients{i},'.mat']);
        for j=1:length(data.p_im)
            temp{j} = [temp{j};...
                data.p_im{j}(data.tight_liver_mask==1)];
        end
    end

    for j=1:length(data.p_im)
       range{j}=prctile(temp{j},[1,99]); 
    end
    clear temp

    
    disp('Computing Haralick texture features (slow)...');
    parfor i=1:length(patients)
        if exist([working_dir,'/norm_data_',patients{i},'.mat'],'file') == 0
            data = load_wrapper([working_dir,'/data_',patients{i},'.mat']);
            f = load_wrapper([working_dir,'/init_features_',patients{i},'.mat']);
            data.glcm = compute_glcm(data.p_im, f.locations, range);
            save_wrapper(data, [working_dir,'/norm_data_',patients{i},'.mat']);
        end
    end
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


function init_features( patients, working_dir )
%INIT_FEATURES Summary of this function goes here
%   Detailed explanation goes here

    % Generating feature cell arrays
    parfor i=1:length(patients)
        if ~exist([working_dir,'/init_features_',patients{i},'.mat'],'file')
            data = load_wrapper([working_dir,'/data_',patients{i},'.mat']);
            features = struct;
            features.locations = find(data.tight_liver_mask);
            features.labels = zeros(length(features.locations),1);

            for c=1:length(features.locations)
                if(data.necro_seg(features.locations(c))==1)
                    features.labels(c)=3;
                elseif(data.tumor_seg(features.locations(c))==1)
                    features.labels(c)=1;
                elseif(data.vasc_seg(features.locations(c))==1)
                    features.labels(c)=2;
                end
            end

            save_wrapper(features, [working_dir,'/init_features_',patients{i},'.mat']);
        end
    end
end
