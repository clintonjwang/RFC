function [] = dicom_display(image,slice)

max_z = length(image(1,1,:));

figure
colormap(jet)

imagesc(image(:,:,slice));
colorbar

while(true)
    waitforbuttonpress
    c = get(gcf,'CurrentCharacter');
    if(strcmp(c,'u'))
        slice = min(slice+1,max_z);
        
        clf
        imagesc(image(:,:,slice));
        colorbar
       
        
        
    elseif(strcmp(c,'d'))
        slice = max(slice-1,1);
        clf
        
        imagesc(image(:,:,slice));
        colorbar
            
    elseif(strcmp(c,'c'))
        close all
        return
    end
    %disp(slice)
end