function printSTORMPoints(stormPoints,fig)
    zmax = max(stormPoints(:,3))+1;
    zmin = min(stormPoints(:,3))-1;
    cmap = colormap('jet');
    figure(fig)
    hold on
    for i = 1:size(stormPoints,1)
       scatter3(stormPoints(i,1),stormPoints(i,2),stormPoints(i,3),3,cmap(ceil((stormPoints(i,3)-zmin)/(zmax-zmin)*64),:),'fill')
    end
end