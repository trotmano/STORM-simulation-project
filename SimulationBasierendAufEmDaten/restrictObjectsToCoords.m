function selectedObjects = restrictObjectsToCoords(objects,xmin,xmax,ymin,ymax,zmin,zmax)
    selectedObjects = {};
    for i = 1:size(objects,2)
        if sum(objects{i}(:,1)<xmin) > 0 || sum(objects{i}(:,1)>xmax)>0 || sum(objects{i}(:,2)<ymin) > 0 || sum(objects{i}(:,2)>ymax)>0 || sum(objects{i}(:,3)<zmin) > 0 || sum(objects{i}(:,3)>zmax)>0
        else
            selectedObjects{1,size(selectedObjects,2)+1} = objects{i};
        end
        
    end


end