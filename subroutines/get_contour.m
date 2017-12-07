%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [contour] = get_contour(mask)% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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