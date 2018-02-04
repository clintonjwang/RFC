function generate_training_data(patients, working_dir)
%GENERATE_TRAINING_DATA Summary of this function goes here

    disp('Extracting training data...');
    f = load_wrapper([working_dir,'/features_',patients{1},'.mat']);
    nf = size(f.intensities,2);
    
    labels = cell(1,length(patients));
    parfor i=1:length(patients)
        labels{i} = load_wrapper([working_dir,'/features_',patients{i},'.mat']);
        labels{i} = labels{i}.labels;
    end
    
    [voxel_data.labels,voxel_data.patients,voxel_data.voxel_locs] = ...
        extract_training_data(labels, 1:length(patients));

    save_wrapper(train, [working_dir,'/voxel_data.mat']);
    
%     voxel_data.features = zeros(length(voxel_data.labels),nf);
%     for i=1:length(patients)
%         loc_i = find(voxel_data.patients==i);
% 
%         if(~isempty(loc_i))
%             f = load_wrapper([working_dir,'/features_',patients{i},'.mat']);
% 
%             for li=loc_i'
%                 ind = voxel_data.voxel_locs(li);
%                 voxel_data.features(li,:) = f.intensities(ind,:);
%             end
%         end
%     end
% 
%     save_wrapper(train, [working_dir,'/voxel_data.mat']);
end

function [labels,patients,locs] = extract_training_data(features,train_indices)
    num_train_pats = length(train_indices); 

    % get median number of voxels in each class
    for i=1:num_train_pats
        normal_temp(i) = length(find(features{train_indices(i)}==0));
        cancer_temp(i) = length(find(features{train_indices(i)}==1));
        vessel_temp(i) = length(find(features{train_indices(i)}==2));
    end
    normal_med = round(median(normal_temp));
    cancer_med = round(median(cancer_temp));
    vessel_med = round(median(vessel_temp)); 
    necrosis_med = 1000;

    % create single list of all voxel indices and patients for each class
    normal_locs=[];
    cancer_locs=[]; 
    vessel_locs=[]; 
    necrosis_locs=[];
    normal_patients=[]; 
    cancer_patients=[]; 
    vessel_patients=[]; 
    necrosis_patients=[];
    for train_indx=1:num_train_pats
        nl = find(features{train_indices(train_indx)}==0);
        cl = find(features{train_indices(train_indx)}==1);
        vl = find(features{train_indices(train_indx)}==2); 
        necl = find(features{train_indices(train_indx)}==3);

        if(length(nl)>normal_med)
            nl = nl(randsample(length(nl),normal_med));
        end
        if(length(cl)>cancer_med)
            cl = cl(randsample(length(cl),cancer_med));
        end
        if(length(vl)>vessel_med)
            vl = vl(randsample(length(vl),vessel_med));
        end
        if(length(necl)>necrosis_med)
            necl = necl(randsample(length(necl),necrosis_med));
        end

        normal_locs = [normal_locs;nl];
        cancer_locs = [cancer_locs;cl];
        vessel_locs = [vessel_locs;vl]; 
        necrosis_locs = [necrosis_locs;necl];

        normal_patients = [normal_patients;train_indices(train_indx)*ones(length(nl),1)];
        cancer_patients = [cancer_patients;train_indices(train_indx)*ones(length(cl),1)];
        vessel_patients = [vessel_patients;train_indices(train_indx)*ones(length(vl),1)];
        necrosis_patients = [necrosis_patients;train_indices(train_indx)*ones(length(necl),1)];
    end
    
    normal_locs = uint32(normal_locs);
    cancer_locs = uint32(cancer_locs);
    vessel_locs = uint32(vessel_locs);
    necrosis_locs = uint32(necrosis_locs);
    normal_patients = uint8(normal_patients);
    cancer_patients = uint8(cancer_patients);
    vessel_patients = uint8(vessel_patients);
    necrosis_patients = uint8(necrosis_patients);

    % T_par=200000; 
    % T_other=50000;
%     T_par=20000;
%     T_other=5000;
% 
%     p = randsample(length(normal_locs),T_par); 
%     normal_locs = normal_locs(p); 
%     normal_patients = normal_patients(p);
% 
%     p = randsample(length(cancer_locs),T_other); 
%     cancer_locs=cancer_locs(p); 
%     cancer_patients = cancer_patients(p);
% 
%     p = randsample(length(vessel_locs),T_other); 
%     vessel_locs=vessel_locs(p); 
%     vessel_patients = vessel_patients(p);
%     labels = [0*ones(T_par,1);1*ones(T_other,1);2*ones(T_other,1);3*ones(length(necrosis_patients),1)];
%     patients = [normal_patients;cancer_patients;vessel_patients;necrosis_patients];
%     locs = [normal_locs;cancer_locs;vessel_locs;necrosis_locs];
    labels = uint8([0*ones(length(normal_locs),1);1*ones(length(cancer_locs),1);2*ones(length(vessel_locs),1);3*ones(length(necrosis_patients),1)]);
    patients = [normal_patients;cancer_patients;vessel_patients;necrosis_patients];
    locs = [normal_locs;cancer_locs;vessel_locs;necrosis_locs];
end