%{
Function: 'flip_image'. 
This function manipulates an image by flipping then transposing it,
slice-by-slice. 

INPUTS: 
    -'image'; the image to manipulate 

OUTPUTS: 
    -'image_o'; the flipped & transposed image 

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function image_o = flip_image(image)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:size(image,3) 
    image_o(:,:,k) = fliplr(fliplr(image(:,:,k))'); 
end

return 