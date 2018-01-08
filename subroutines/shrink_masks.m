function data = shrink_masks(data, N1, N2, N3, train_bool)
%SHRINK_MASKS Summary of this function goes here
%   Detailed explanation goes here   
    temp = zeros(size(data.art)); 
    x_max=min(size(data.liver_mask,1),N1); 
    y_max=min(size(data.liver_mask,2),N2);
    z_max=min(size(data.liver_mask,3),N3);
    temp(1:x_max,1:y_max,1:z_max)=data.liver_mask(1:x_max,1:y_max,1:z_max);
    data.iso_full_size = [x_max y_max z_max];
    data.liver_mask=temp;

    if train_bool
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
    end

return