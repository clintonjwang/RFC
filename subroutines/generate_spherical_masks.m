function M = generate_spherical_masks(sl)

A=zeros((2*sl)+1,(2*sl)+1,(2*sl)+1);
center=[sl+1,sl+1,sl+1]; 

for i=1:size(A,1)
    for j=1:size(A,2)
        for k=1:size(A,3)
            A(i,j,k) = norm([i,j,k] - center);  
        end
    end
end

A=round(A); 

ind=1; 
%{
for r=1:sl
    B=(A==r);
    M{ind} = find(B);
    [I,J,K]=ind2sub(size(A),M{ind});
    P=[I,J,K]-repmat(center,numel(I),1);
    
    ind=ind+1;
end
%}

for r=1:sl
    M{ind} = find(A==r);
    ind=ind+1; 
end

%{
for r=1:sl-1
    M{ind} = find(A<=r);
    ind=ind+1; 
end
%}

return 
end