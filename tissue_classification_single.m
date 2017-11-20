function [vasc_mask, nec_mask, viatumor_mask, paren_mask] = ...
        tissue_classification_single( data, f, intensities, train_mode )
%TISSUE_CLASSIFICATION_SINGLE Classify a liver based on...
%   Detailed explanation goes here
dir_main = '.';
dir_scratch60 = 'features';
dir_feat = 'features';
model_dir = 'models';

%initialize global variables
ntrees=800; %trees in each random forest 
RAC=3; %number of iterations of classification
sl=8; %structured context patch size
sl_spherical=5; %spherical context patch size
num_bins=5; %number of histogram bins for spherical histogram context features 
sa=6; %number of auto-context features to sample per dimension
min_leaf_size=50; %random forest parameter: mininum leaf size
num_classes=4; %number of tissue classes
M = generate_spherical_masks(sl_spherical); %generate spherical masks
num_context_features = (sa^3) * (num_classes - 1)... %structured context features
         + ((num_classes-1)*num_bins*length(M)); %plus spherical context features 

%loop through rounds of auto-context
for r=1:RAC
    disp(['Current round of auto-context...',num2str(r)]);

    %clear trees from previous round of label-context
    clear C

%     fileID = fopen([model_dir,'/tree_',num2str(r),'.mat'],'r');
%     C = fread(fileID,[numel(f.labels), num_context_features],'double');
%     fclose(fileID);
    C = load([model_dir,'/tree_',num2str(r),'.mat']);

    %read in appearance feature binary file
%     fileID = fopen([dir_scratch60, '/intensities_',num2str(td),'.bin'],'r');
%     intensities = fread(fileID,[numel(f.labels),...
%         f.num_intensity_features],'double');
%     fclose(fileID);

    %variable holding tissue probability estimates 
    classification{td}=zeros(numel(f.labels),num_classes);

    tic
    %if first round of classification: appearance features only 
    if(r==1)
        for t=1:ntrees
            [~,scores] = predict(C{t},intensities);
            classification{td} = classification{td} + scores;
        end
    %for subsequent rounds of classification, incorporate auto-context 
    %features 
    else
        %read in binary file holding auto-context features for
        fileID = fopen([dir_scratch60, '/context_',num2str(td),'.bin'],'r');
        auto_context_features = fread(fileID,[numel(f.labels),...
            num_context_features],'double');
        fclose(fileID);
        for t=1:ntrees
            [~,scores] = predict(C{t},[intensities,auto_context_features]);
            classification{td} = classification{td} + scores;
        end
    end
    classification{td}=classification{td}./ntrees;
    toc

    %clear the previously trained random forest 
    clear RF

    %compute the auto-context features, unless it is the final iteration
    if(r~=RAC)
        auto_context_features = extract_autocontext_features(...
            classification,data,locations,sl,sa,train.data{cv_iteration},...
            train.train_patients{cv_iteration},train.train_locs{cv_iteration},train_indices{cv_iteration},...
            num_context_features,sl_spherical,num_bins,dir_scratch60);
    else
        [vasc_mask, nec_mask, viatumor_mask, paren_mask] = get_masks(classification)
        write_ids_mask(vasc_mask, '.', '.', 'vasculature_mask');
        write_ids_mask(nec_mask, '.', '.', 'necrosis_mask');
        write_ids_mask(viatumor_mask, '.', '.', 'viable_tumor_mask');
        write_ids_mask(paren_mask, '.', '.', 'parenchyma_mask');
    end
end

return
end




%dummy function to allow binary file saves inside parfor-loop 
function [] = save_dummy_binary(A,dir_scratch60)

fileID = fopen([dir_scratch60, '/context.bin'],'w');
fwrite(fileID,A,'double');
fclose(fileID);

return
end


%function to compute auto-context features 
function auto_context_features = extract_autocontext_features(classification,...
    data,locations,sl,sa,train_data,train_patients,...
    train_locs,train_indices,num_context_features,sl_spherical,num_bins,save_dir,...
    f)

%initialize variables
num_classes = 4;
num_train_points = length(train_patients);

disp('Extracting auto-context features...')

T=[];
acf=[];
    
%compute label context features
tic
maps = generate_maps_dummy(classification,f.sz,locations,num_classes);

auto_context_features = [compute_patches(maps,locations,sl,data,sa),...
    compute_patches_spherical(maps,locations,sl_spherical,data,num_bins)];
toc

%save auto-context features in binary file
save_dummy_binary(auto_context_features,save_dir);

%augment the appearance features in the training data with auto-context
%features
temp=zeros(num_train_points,num_context_features);
for td=train_indices
    temp(T{td},:) = acf{td};
end

train_data_ac = [train_data,temp];

return
end


%compute structured auto-context features 
function temp = compute_patches(maps,locations,sl,data_i,sa)
%initialize variables...
sz=size(data_i.tight_liver_mask);
[ii,jj,kk] = ind2sub(sz,locations);

for c=1:size(maps,1)
    for c1=1:size(maps,2)
        patches{c,c1}=zeros(length(ii),sa(c1)^3);
    end
end
%loop through data points
for d=1:length(ii)
    for c1=1:size(maps,2)
        
        x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,sz(1)),sa(c1)));
        y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,sz(2)),sa(c1)));
        z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,sz(3)),sa(c1)));
        sample_locs = combvec(x,y,z);
        
        inds = sub2ind(sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));
        
        for c=1:size(maps,1)
            patches{c,c1}(d,1:(sa(c1))^3) = maps{c,c1}(inds);
        end
   
    end
end

temp=[];
for c=1:size(maps,1)
    for c1=1:size(maps,2)
        temp=[temp,patches{c,c1}];
    end
end

return
end


%compute spherical histogram auto-context features 
function temp = compute_patches_spherical(maps,locations,sl,data_i,num_bins)

M = generate_spherical_masks1(sl);
num_spheres=length(M);
sz=size(data_i.tight_liver_mask); 
[ii,jj,kk] = ind2sub(sz,locations);
num_maps = length(maps);
N=(2*sl)+1;

edges = [linspace(0,1,num_bins+1)];

for c=1:length(maps)
    patches{c}=zeros(length(ii),(num_spheres*num_bins)); 
end

for d=1:length(ii)
    
    x=round(linspace(max(ii(d)-sl,1),min(ii(d)+sl,sz(1)),N));
    y=round(linspace(max(jj(d)-sl,1),min(jj(d)+sl,sz(2)),N));
    z=round(linspace(max(kk(d)-sl,1),min(kk(d)+sl,sz(3)),N));
    sample_locs = combvec(x,y,z);

    inds = sub2ind(sz,sample_locs(1,:),sample_locs(2,:),sample_locs(3,:));
    
    for c=1:length(maps)
        for s=1:num_spheres
            h=histcounts(maps{c}(inds(M{s})),edges);
            patches{c}(d,((s-1)*num_bins)+1:s*num_bins) = h;
        end
    end
    
end

temp=[];
for c=1:num_maps
   temp=[temp,patches{c}];  
end

return
end
