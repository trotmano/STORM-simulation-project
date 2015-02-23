function dists = getDistanceToFluorophor(posFluorophor, idx, idxF)
    dists = [];
    for i = 1:size(posFluorophor,1)
        minDistRight = 1e9;
        minDistWrong = 1e9;
        for j = 1: size(posFluorophor,1)
            tmpDist = norm(posFluorophor(j,:)- posFluorophor(i,:));
            if idx(i) == idx(j) && tmpDist<minDistRight && tmpDist>0
                minDistRight = tmpDist;
            elseif idx(i)~=idx(j) && tmpDist<minDistWrong
                minDistWrong = tmpDist;
            end
            
            
        end
        dists(idxF==idx(i),1) = minDistRight;
        dists(idxF==idx(i),2) = minDistWrong;
    end
end

