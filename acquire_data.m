% textscan(fopen('Test.txt','r'),'%s','Delimiter','\n')

%input 'patients' should be cell array of patient IDs
function data = acquire_data(patients)

disp('Acquiring patient data...');
dest = './training_data/'

if exist([dest,patients{1},'/nii_files/whole_liver.nii'],'file')>0
    disp('Data already acquired!')
    return
end

%parfor i=1:length(patients)
for i=1:length(patients)
    data{i} = get_liver_data(patients{i}, dest);
    disp([patients{i}, ' acquired!']);
end

return

function data = get_liver_data(patient, dest)

%directory containing folders with patient data 
%dest='C:/Users/John/Desktop/Data/TissueClassification/';

    
disp(patient);
data.patID = patient;

R=1.75; %desired low-resolution in mm

%reslice the pre, 20s, 70s, to make them isotropic
temp = load_nii([dest,patient,'/nii_files/20s.nii.gz']);
temp_res(1) = temp.hdr.dime.pixdim(2); %extract pixel volume
temp_res(2) = temp.hdr.dime.pixdim(3); % dido
temp_res(3) = temp.hdr.dime.pixdim(4); % dido

[N1,N2,N3] = size(double(flip_image(temp.img)));

%make nii file from whole liver segmentation 
fpath = [dest,patient,'/segs/whole_liver.ids'];
liver_mask = get_mask(fpath,N1,N2,N3);
liver_nii = make_nii(flip_image(liver_mask),temp_res);
save_nii(liver_nii,[dest,patient,'/nii_files/whole_liver.nii']);

%flags indicating the existence of vessel, tumor, and necrosis
%segmentations 
vessel_exist = exist([dest,patient,'/segs/vessel.ids'],'file')>0;
tumor_exist = exist([dest,patient,'/segs/tumor.ids'],'file')>0;
necrosis_exist = exist([dest,patient,'/segs/nec.ids'],'file')>0;

%make nii file from vessel segmentation 
if(vessel_exist)
    fpath = [dest,patient,'/segs/vessel.ids'];
    vessel_mask = get_mask(fpath,N1,N2,N3);
else
    vessel_mask = zeros(N1,N2,N3);
end
vessel_nii = make_nii(flip_image(vessel_mask),temp_res);
save_nii(vessel_nii,[dest,patient,'/nii_files/vessel.nii']);

%make nii file from tumor segmentation 
if(tumor_exist)
    fpath = [dest,patient,'/segs/tumor.ids'];
    tumor_mask = get_mask(fpath,N1,N2,N3);
else
    tumor_mask = zeros(N1,N2,N3);
end
tumor_nii = make_nii(flip_image(tumor_mask),temp_res);
save_nii(tumor_nii,[dest,patient,'/nii_files/tumor.nii']);

%make nii file from tumor segmentation 
if(necrosis_exist)
    necrosis_mask = get_mask([dest,patient,'/segs/nec.ids'],N1,N2,N3);
else
    necrosis_mask = zeros(N1,N2,N3);
end
necrosis_nii = make_nii(flip_image(necrosis_mask),temp_res);
save_nii(necrosis_nii,[dest,patient,'/nii_files/necrosis.nii']);

%reslice 20s image to be isotropic 
reslice_nii([dest,patient,'/nii_files/20s.nii.gz'],...
    [dest,patient,'/nii_files/20s_isotropic.nii.gz'],...
    [R,R,R]);

%reslice pre image to be isotropic
reslice_nii([dest,patient,'/nii_files/pre_reg.nii'],...
    [dest,patient,'/nii_files/pre_reg_isotropic.nii'],...
    [R,R,R]);

%reslice 70s to be isotropic 
reslice_nii([dest,patient,'/nii_files/70s_reg.nii'],...
    [dest,patient,'/nii_files/70s_reg_isotropic.nii'],...
    [R,R,R]);

%reslice whole liver segmentation to be isotropic 
reslice_nii([dest,patient,'/nii_files/whole_liver.nii'],...
    [dest,patient,'/nii_files/whole_liver_isotropic.nii'],...
    [R,R,R]);

%reslice vessel segmentation to be isotropic 
reslice_nii([dest,patient,'/nii_files/vessel.nii'],...
    [dest,patient,'/nii_files/vessel_isotropic.nii'],...
    [R,R,R]);

%reslice necrosis segmentation to be isotropic 
reslice_nii([dest,patient,'/nii_files/necrosis.nii'],...
    [dest,patient,'/nii_files/necrosis_isotropic.nii'],...
    [R,R,R]);

%reslice tumor segmentation to be isotropic 
reslice_nii([dest,patient,'/nii_files/tumor.nii'],...
    [dest,patient,'/nii_files/tumor_isotropic.nii'],...
    [R,R,R]);

%reslice t2 image bias field estimate to be isotropic 
reslice_nii([dest,patient,'/nii_files/t2_bfc_reg.nii'],...
    [dest,patient,'/nii_files/t2_bfc_reg_isotropic.nii'],...
    [R,R,R]);

%load pre image 
data.pre = load_nii([dest,patient,'/nii_files/pre_reg_isotropic.nii']);
data.pre = double(flip_image(data.pre.img));

%load 20s image 
data.art = load_nii([dest,patient,'/nii_files/20s_isotropic.nii.gz']);
data.resolution(1) = data.art.hdr.dime.pixdim(2); %extract pixel volume
data.resolution(2) = data.art.hdr.dime.pixdim(3); % dido
data.resolution(3) = data.art.hdr.dime.pixdim(4); % dido
data.art = double(flip_image(data.art.img));

%load 70s image 
data.pv = load_nii([dest,patient,'/nii_files/70s_reg_isotropic.nii']);
data.pv = double(flip_image(data.pv.img));

%load t2 image
data.t2 = load_nii([dest,patient,'/nii_files/t2_bfc_reg_isotropic.nii']);
data.t2 = double(flip_image(data.t2.img));

%load T1 bias field estimate 
data.bf = load_nii([dest,patient,'/nii_files/bias_field_isotropic.nii']);
data.bf = double(flip_image(data.bf.img));

%get T1 image dimensions
[N1,N2,N3] = size(data.art);

%load liver_mask
data.liver_mask = load_nii([dest,patient,'/nii_files/whole_liver_isotropic.nii']);
data.liver_mask = double(flip_image(data.liver_mask.img));
data.liver_mask = data.liver_mask>0;

%load vessel mask 
data.vessel_mask = load_nii([dest,patient,'/nii_files/vessel_isotropic.nii']);
data.vessel_mask = double(flip_image(data.vessel_mask.img));
data.vessel_mask = data.vessel_mask>0;

%load necrosis mask 
data.necrosis_mask = load_nii([dest,patient,'/nii_files/necrosis_isotropic.nii']);
data.necrosis_mask = double(flip_image(data.necrosis_mask.img));
data.necrosis_mask = data.necrosis_mask>0;

%load tumor mask 
data.tumor_mask = load_nii([dest,patient,'/nii_files/tumor_isotropic.nii']);
data.tumor_mask = double(flip_image(data.tumor_mask.img));
data.tumor_mask = data.tumor_mask>0;



%shrink the masks to proper size if necessary
temp = zeros(size(data.art)); 
x_max=min(size(data.liver_mask,1),N1); 
y_max=min(size(data.liver_mask,2),N2);
z_max=min(size(data.liver_mask,3),N3);
temp(1:x_max,1:y_max,1:z_max)=data.liver_mask(1:x_max,1:y_max,1:z_max);  
data.liver_mask=temp; 

temp = zeros(size(data.art)); 
x_max=min(size(data.vessel_mask,1),N1); 
y_max=min(size(data.vessel_mask,2),N2);
z_max=min(size(data.vessel_mask,3),N3);
temp(1:x_max,1:y_max,1:z_max)=data.vessel_mask(1:x_max,1:y_max,1:z_max);  
data.vessel_mask=temp;

temp = zeros(size(data.art)); 
x_max=min(size(data.necrosis_mask,1),N1); 
y_max=min(size(data.necrosis_mask,2),N2);
z_max=min(size(data.necrosis_mask,3),N3);
temp(1:x_max,1:y_max,1:z_max)=data.necrosis_mask(1:x_max,1:y_max,1:z_max);  
data.necrosis_mask=temp;

temp = zeros(size(data.art)); 
x_max=min(size(data.tumor_mask,1),N1); 
y_max=min(size(data.tumor_mask,2),N2);
z_max=min(size(data.tumor_mask,3),N3);
temp(1:x_max,1:y_max,1:z_max)=data.tumor_mask(1:x_max,1:y_max,1:z_max);  
data.tumor_mask=temp;

%compute tightened liver mask
r=4;  
rs=r/data.resolution(1); 
D=bwdist(1-data.liver_mask);
data.tight_liver_mask = zeros(N1,N2,N3);
data.tight_liver_mask(D>rs) = data.liver_mask(D>rs);

data.necrosis_mask = data.necrosis_mask .* data.tight_liver_mask;
data.tumor_mask = data.tumor_mask .* data.tight_liver_mask .*(1-data.necrosis_mask);
data.vessel_mask = data.vessel_mask.*data.tight_liver_mask.*(1-data.tumor_mask).*(1-data.necrosis_mask);

%find edges of the liver and crop the images to this dimension
[i,j,k] = ind2sub(size(data.liver_mask),find(data.liver_mask));
i_min = min(i);
i_max = max(i);

j_min=min(j);
j_max=max(j);

k_min=min(k);
k_max=max(k);

sf={'pre','art','pv','t2','liver_mask','tumor_mask','bf',...
    'tight_liver_mask','vessel_mask','necrosis_mask'};

for sf_count=1:length(sf)
    if(isfield(data,sf{sf_count}))
        data.(sf{sf_count}) = data.(sf{sf_count})(i_min:i_max, ...
            j_min:j_max,k_min:k_max);
    end
end

%get contours for later visualization
data.liver_contour = get_contour(data.liver_mask);
data.tight_liver_contour = get_contour(data.tight_liver_mask);
data.tumor_contour = get_contour(data.tumor_mask);
data.vessel_contour = get_contour(data.vessel_mask);
data.necrosis_contour = get_contour(data.necrosis_mask);

return

function [mask] = get_mask(fpath,N1,N2,N3)%

fileID = fopen(fpath);
A = fread(fileID);
ind = find(A);
[i, j, k]=ind2sub([N2,N1,N3],ind);
fclose('all');

mask = zeros(N1,N2,N3);
for count=1:length(i)
    mask(i(count),j(count),k(count))=1;
end

mask = transpose_mask_slices(mask);

return

function [contour] = get_contour(mask)

[N1,N2,N3] = size(mask);

for n=1:N3
    perim(:,:,n) = bwperim(mask(:,:,n));
end

for n3=1:N3
    count=1;
    coordinates=[];
    for n1=1:N1
        for n2=1:N2
            if(perim(n1,n2,n3)==1)
                coordinates(count,1)=n1;
                coordinates(count,2)=n2;
                count = count + 1;
            end
        end
    end
    contour{n3} = coordinates;
end

return

function mask1 = transpose_mask_slices(mask)

[N1,N2,N3] = size(mask);

mask1 = zeros(N1,N2,N3);

for k=1:N3
    for i=1:N1
        for j=1:N2
            if(mask(i,j,k)==1)
                mask1(j,i,k)=1;
            end
        end
    end
end

return

function image_o = flip_image(image)%

for k=1:size(image,3)
    image_o(:,:,k) = fliplr(fliplr(image(:,:,k))');
end

return