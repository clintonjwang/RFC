function [mask] = get_mask(fpath)

fileID = fopen([fpath,'.ics']);
A = fscanf(fileID,'%c');
D = sscanf(A(regexp(A, 'sizes\t')+8:regexp(A, 'layout\tcoordinates')-3), '%d');
N1 = D(1);
N2 = D(2);
N3 = D(3);

fileID = fopen([fpath,'.ids']);
A = fread(fileID);
ind = find(A);
[i, j, k]=ind2sub([N2,N1,N3],ind);
fclose('all'); 

mask = zeros(N1,N2,N3);
for count=1:length(i) 
    mask(i(count),j(count),k(count))=1; 
end

mask = transpose_mask_slices(mask, 'r');

return 