function [] = liver_display(image,pat,slice,tt)

%image = image .* pat.tight_liver_mask;


liver_contour = pat.tight_liver_contour;
if(tt==1)
    cancer_contour = pat.tumor_contour;
    vessel_contour = pat.vessel_contour;
    necrosis_contour = pat.necrosis_contour;
end

max_z = length(image(1,1,:));

figure

colormap(gray)

imagesc(image(:,:,slice));

hold on

if(tt==1)
    if(~isempty(necrosis_contour{slice}))
        disp('flag')
        scatter(necrosis_contour{slice}(:,2),necrosis_contour{slice}(:,1),5,'b');
    end
    
    if(~isempty(cancer_contour{slice}))
        scatter(cancer_contour{slice}(:,2),cancer_contour{slice}(:,1),5,'r');
    end
    
    if(~isempty(vessel_contour{slice}))
        scatter(vessel_contour{slice}(:,2),vessel_contour{slice}(:,1),5,'y');
    end
end

if(~isempty(liver_contour{slice}))
    scatter(liver_contour{slice}(:,2),liver_contour{slice}(:,1),5,'g');
end

colorbar

while(true)
    waitforbuttonpress
    c = get(gcf,'CurrentCharacter');
    if(strcmp(c,'u'))
        slice = min(slice+1,max_z);
        
        clf
        imagesc(image(:,:,slice));
        
        hold on
        if(tt==1)
            if(~isempty(necrosis_contour{slice}))
                scatter(necrosis_contour{slice}(:,2),necrosis_contour{slice}(:,1),5,'b');
            end
            
            if(~isempty(cancer_contour{slice}))
                scatter(cancer_contour{slice}(:,2),cancer_contour{slice}(:,1),5,'r');
            end
            
            
            if(~isempty(vessel_contour{slice}))
                scatter(vessel_contour{slice}(:,2),vessel_contour{slice}(:,1),5,'y');
            end
        end
        
        
        if(~isempty(liver_contour{slice}))
            scatter(liver_contour{slice}(:,2),liver_contour{slice}(:,1),5,'g');
        end
        
        colorbar
        
    elseif(strcmp(c,'d'))
        slice = max(slice-1,1);
        clf
        
        imagesc(image(:,:,slice));
        
        hold on
        
        if(tt==1)
            if(~isempty(necrosis_contour{slice}))
                scatter(necrosis_contour{slice}(:,2),necrosis_contour{slice}(:,1),5,'b');
            end
            
            if(~isempty(cancer_contour{slice}))
                scatter(cancer_contour{slice}(:,2),cancer_contour{slice}(:,1),5,'r');
            end
            
            
            
            if(~isempty(vessel_contour{slice}))
                scatter(vessel_contour{slice}(:,2),vessel_contour{slice}(:,1),5,'y');
            end
            
        end
        
        
        if(~isempty(liver_contour{slice}))
            scatter(liver_contour{slice}(:,2),liver_contour{slice}(:,1),5,'g');
        end
        
        
        colorbar
    elseif(strcmp(c,'c'))
        close all
        return
    end
    %disp(slice)
end