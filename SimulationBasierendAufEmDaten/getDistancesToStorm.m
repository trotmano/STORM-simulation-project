function dists = getDistancesToStorm(objects,stormPoints, idxF)
    dists = [];
    for i = 1:size(stormPoints,1)
        for j = 1: size(stormPoints,1)
            minDistRight = 1e9;
            minDistWrong = 1e9;
            for j = 1:size(objects,2)
                tmpDist = norm(stormPoints(j,:)- stormPoints(i,:));
                if idxF(i) == idxF(j) && tmpDist<minDistRight
                    minDistRight = tmpDist;
                elseif idxF(i)~=idxF(j) && tmpDist<minDistWrong
                    minDistWrong = tmpDist;
                end
            end
            dists(i,:) = [minDistRight, minDistWrong];
        end
        
    end
end

