function maps = generate_maps_dummy(class,sz,locs,num_classes)

%sigma=[1,2];

for c=2:num_classes
    maps{c-1,1} = generate_classification_visualization(class(:,c),sz,locs);
    %{
    for j=1:numel(sigma)
        maps{c-1,j+1} = imgaussfilt3(maps{c-1,1},sigma(j));
    end
    %}
end
return
end

function map = generate_classification_visualization(p,sz,locs) 
 
map=zeros*ones(sz);

for i=1:length(locs)
    map(locs(i))=p(i); 
end

map(1,1,:) = ones(1,size(map,3));

end