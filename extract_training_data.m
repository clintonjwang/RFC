function [train_labels,train_patients,train_locs] = ...
    extract_training_data(features,train_indices)
% Does not modify features

num_train_pats = length(train_indices); 

for i=1:num_train_pats
    normal_temp(i) = length(find(features{train_indices(i)}==0));
    cancer_temp(i) = length(find(features{train_indices(i)}==1));
    vessel_temp(i) = length(find(features{train_indices(i)}==2));
end

normal_med = round(median(normal_temp));
cancer_med = round(median(cancer_temp));
vessel_med = round(median(vessel_temp)); 
necrosis_med = 1000;

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

% T_par=200000; 
% T_other=50000;
T_par=20000; 
T_other=5000;

p=randsample(length(normal_locs),T_par); 
normal_locs=normal_locs(p); 
normal_patients = normal_patients(p);

p=randsample(length(cancer_locs),T_other); 
cancer_locs=cancer_locs(p); 
cancer_patients = cancer_patients(p);

p=randsample(length(vessel_locs),T_other); 
vessel_locs=vessel_locs(p); 
vessel_patients = vessel_patients(p);

train_labels = [0*ones(T_par,1);1*ones(T_other,1);2*ones(T_other,1);3*ones(length(necrosis_patients),1)];
train_patients = [normal_patients;cancer_patients;vessel_patients;necrosis_patients];
train_locs = [normal_locs;cancer_locs;vessel_locs;necrosis_locs]; 

return
end