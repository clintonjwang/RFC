function map = generate_classification_visualization(p,sz,locs) 
 
map=zeros*ones(sz); 

for i=1:length(locs)
    map(locs(i))=p(i); 
end

map(1,1,:) = ones(1,size(map,3));

end