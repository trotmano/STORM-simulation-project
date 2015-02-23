function printSTORMPointsClassified(stormPoints, dists, distsStorm, fig)
    zmax = max(stormPoints(:,3))+1;
    zmin = min(stormPoints(:,3))-1;
    cmap = colormap('jet');
    figure(fig)
    hold on
    factor = 5;
    pointSize = 10;
    for i = 1:size(stormPoints,1)
        if factor*dists(i,1)<dists(i,2)
            col1 = 0;
        elseif dists(i,1)>dists(i,2)
            col1 = 1;
        else 
            col1 = (dists(i,1)*factor - dists(i,2))/(factor-1)/dists(i,2);
        end
        
        if factor*distsStorm(i,1)<distsStorm(i,2)
            col2 = 0;
        elseif distsStorm(i,1)>distsStorm(i,2)
            col2 = 1;
        else 
            col2 = (distsStorm(i,1)*factor - distsStorm(i,2))/(factor-1)/distsStorm(i,2);
        end
        
        
        color = [(col1+col2),1-(col1+col2),0];
        scatter3(stormPoints(i,1),stormPoints(i,2),stormPoints(i,3),pointSize,color,'fill')
    end
end