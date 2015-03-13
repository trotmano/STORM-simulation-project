function [listStartPoints, listEndPoints,idxs] = findLines(points, bspnm, pabs,aoa,loa,rof)
    nbrPoints = size(points,2);
    listStartPoints = [];
    listEndPoints = [];
    idxs = []; %idx contains information from which structure the Antibody originates
    for i = 1:nbrPoints
        [lengthOfStructure, cummulativeLengths] = getLengthOfStructure(points{1,i});
        aoa = randn(floor(bspnm*lengthOfStructure),1)*20/180*pi+(90/180*pi);
        for j = 1:floor(bspnm*lengthOfStructure) %length of line
            
            randomNumber = rand(1,1);
            if randomNumber < pabs
                
                idx = sum(cummulativeLengths<(j/bspnm))+1;
                lineVec = [ points{1,i}(idx+1,1) - points{1,i}(idx,1),...
                            points{1,i}(idx+1,2)-points{1,i}(idx,2),...
                            points{1,i}(idx+1,3)-points{1,i}(idx,3)];
                alpha = rand(1)*2*pi;
                if (size(aoa,1)>1)
                    vecOrth = getVector(aoa(j), rof,alpha);
                    vec = getVector(aoa(j), loa,alpha);
                else
                    vecOrth = getVector(aoa, rof,alpha);
                    vec = getVector(aoa, loa,alpha);
                end
                rotVec = findRotation(vec,points{1,i}(idx:idx+1,:));
                rotVecOrth = findRotation(vecOrth, points{1,i}(idx:idx+1,:));
                startPoint = points{1,i}(idx+1,:)+((j-1)/bspnm-cummulativeLengths(idx))*lineVec/norm(lineVec)+rotVecOrth';
                endPoint = points{1,i}(idx+1,:)+((j-1)/bspnm-cummulativeLengths(idx))*lineVec/norm(lineVec)+rotVecOrth'+rotVec';
                if isnan(rotVec(1))
                else
                    listStartPoints(size(listStartPoints,1)+1,:) = startPoint;
                    listEndPoints(size(listEndPoints,1)+1,:) = endPoint;
                    idxs(size(idxs,1)+1,:) = i;
                end
            end
        end
    end    
end