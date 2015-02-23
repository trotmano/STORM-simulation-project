function dists = getDistancesToStructures(objects,stormPoints, idxF)
    dists = [];
    for i = 1:size(stormPoints,1)
        minDistRight = 1e9;
        minDistWrong = 1e9;
        for j = 1:size(objects,2)
            tmpDist = getDistanceToStructure(objects{j}, stormPoints(i,:));
            if idxF(i) == j && tmpDist<minDistRight
                minDistRight = tmpDist;
            elseif idxF(i)~=j && tmpDist<minDistWrong
                minDistWrong = tmpDist;
            end
        end
        dists(i,:) = [minDistRight, minDistWrong];
    end
end

function dist = getDistanceToStructure(object, stormPoint)
    dist = 1e9; %http://geomalgorithms.com/a02-_lines.html
    w = stormPoint-object(1,:);
    v = object(2,:) - object(1,:);
    b = w.*v/(v.*v);
    pb = object(1,:) + b * v;
    if b<0
        dist = norm(stormPoint- object(1,:));
    elseif b > norm(v)
        dist = norm(stormPoint- object(2,:));
    else
        dist = norm(stormPoint- pb);
    end
         
end