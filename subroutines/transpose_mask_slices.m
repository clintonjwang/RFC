function out_mask = transpose_mask_slices(mask, mode)
%Required to read/write .ids files

if mode == 'r'
    [N1,N2,N3] = size(mask);

    out_mask = zeros(N1,N2,N3);

    for k=1:N3
        for i=1:N1
            for j=1:N2
                if(mask(i,j,k)==1)
                    out_mask(j,i,k)=1;
                end
            end
        end
    end
    
elseif mode == 'w'
    [N2,N1,N3] = size(mask);

    out_mask = zeros(N1,N2,N3);

    for k=1:N3
        for j=1:N2
            for i=1:N1
                if(mask(j,i,k)==1)
                    out_mask(i,j,k)=255;
                end
            end
        end
    end
    
else
    disp('Invalid mode in transpose_mask_slices')
end
return