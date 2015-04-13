function [stormPoints, idxF, idxST] = findStormPoints(listEndPoints, abpf, sxy, sz, bd,fpab,background)
    stormPoints = [];
    idxF = []; %idxF contains the information to which structure the fluorophore belongs
    idxST = []; %idxST contains the information to which structure each localization belongs
%     for i = 1:size(listEndPoints,1)
%             for j = 1:abpf
%                 x= listEndPoints(i,1)+randn(1,1) * sxy;
%                 y= listEndPoints(i,2)+randn(1,1) * sxy;
%                 z= listEndPoints(i,3)+randn(1,1) * sz;
%                 stormPoints(size(stormPoints,1)+1,1) = x;
%                 stormPoints(size(stormPoints,1),2) = y;
%                 stormPoints(size(stormPoints,1),3) = z;
%                 idxST(size(idxST,1)+1,:) = idx(i);
%                 idxF(size(idxF,1)+1,:) = i;
%             end
%     end
    if background %unspecific labeling
        ilpmm3 = 50; %incorrect localizations per micrometer ^3
        xmin = min(listEndPoints(:,1));
        xmax = max(listEndPoints(:,1));
        ymin = min(listEndPoints(:,2));
        ymax = max(listEndPoints(:,2));
        zmin = min(listEndPoints(:,3));
        zmax = max(listEndPoints(:,3));
        numberOfIncorrectLocalizations = floor(ilpmm3*(xmax-xmin)/1e3*(ymax-ymin)/1e3*(zmax-zmin)/1e3);
        x = rand(numberOfIncorrectLocalizations,1) * (xmax -xmin) + xmin;
        y = rand(numberOfIncorrectLocalizations,1) * (ymax -ymin) + ymin;
        z = rand(numberOfIncorrectLocalizations,1) * (zmax -zmin) + zmin;
        listEndPoints = [listEndPoints;[x,y,z]];
    end
    
    %add additional fluorophores near the endpoint
    if fpab~= 1
        idx = abs(floor(randn(size(listEndPoints,1),1)*fpab+fpab));
        idx(idx==0) = 1;
        listEndPointsAugmented=[];
        for i=1:max(idx)
            alteredPoints = listEndPoints(idx>=i,:);
            alteredPoints = alteredPoints+randn(size(alteredPoints))*3;
            listEndPointsAugmented = [listEndPointsAugmented;alteredPoints];
        end
        listEndPoints = listEndPointsAugmented;
    end    
    
    
    
    nbrBlinkingEvents = randn(size(listEndPoints,1),1)*sqrt(abpf) + abpf;
    nbrBlinkingEvents(nbrBlinkingEvents<0) = 0;
    for i = 1:floor(max(max(nbrBlinkingEvents)))
        idx = nbrBlinkingEvents >= i;
        x = listEndPoints(idx,1) + randn(size(listEndPoints(idx),1),1) * sxy;
        y = listEndPoints(idx,2) + randn(size(listEndPoints(idx),1),1) * sxy;
        z = listEndPoints(idx,3) + randn(size(listEndPoints(idx),1),1) * sz;
        intensities = abs(random('exp',0.0065,[size(listEndPoints(idx),1),1]))*3.8e4+1e3;
        stormPoints = [stormPoints;[x,y,z, intensities]];
    end
    
%     for i = 1:abpf
%         x = listEndPoints(:,1) + randn(size(listEndPoints,1),1) * sxy;
%         y = listEndPoints(:,2) + randn(size(listEndPoints,1),1) * sxy;
%         z = listEndPoints(:,3) + randn(size(listEndPoints,1),1) * sz;
%         intensity = abs(randn(size(listEndPoints,1),1))*3.8e4;
%         stormPoints = [stormPoints;[x,y,z, intensities]];
%         idxST = [idxST;idx];
%         idxF = [idxF,linspace(1,size(listEndPoints,1),size(listEndPoints,1))];
%     end
    if size(stormPoints,1) == 0
        
    else
        fluorophoresPerFrame = (max(stormPoints(:,1))-min(stormPoints(:,1)))...
                            *(max(stormPoints(:,2))-min(stormPoints(:,2)))...
                            *bd;
        if fluorophoresPerFrame<1
          fluorophoresPerFrame = 1;
        end
        stormPoints = [stormPoints,randi([0,ceil(size(stormPoints,1)/fluorophoresPerFrame)],size(stormPoints,1),1)];
        mergedPSFs = 1;
        psfwidth = 200;
        affectingFactor = 2;
        if mergedPSFs
            disp('mergedPSF')
            for i = 1:max(stormPoints(:,5))
                i;
                idx = stormPoints(:,5)==i;
                [idx2,tmp] = find(idx==1);
                dists = pdist2(stormPoints(idx,1:2),stormPoints(idx,1:2)); %only x and y
                dists = dists + tril(ones(size(dists))*9e9,1);
                [firstLoc,secondLoc] = find(dists<psfwidth);
                [firstLoc2, secondLoc2] = find(logical(dists<affectingFactor*psfwidth).*logical(dists>psfwidth));
                meanCoords = (stormPoints(idx2(firstLoc),1:5)+stormPoints(idx2(secondLoc),1:5))/2.;

                diffVec = stormPoints(idx2(firstLoc2),1:2)-stormPoints(idx2(secondLoc2),1:2);
                for j =1:size(diffVec,1)
                    stormPoints(idx2(firstLoc2(j)),1:2)=stormPoints(idx2(firstLoc2(j)),1:2)-(affectingFactor*psfwidth-norm(diffVec(j,:)))/(affectingFactor*psfwidth)*0.5*diffVec(j,:);
                    stormPoints(idx2(secondLoc2(j)),1:2)=stormPoints(idx2(secondLoc2(j)),1:2)+(affectingFactor*psfwidth-norm(diffVec(j,:)))/(affectingFactor*psfwidth)*0.5*diffVec(j,:);
                    lengthDiffVec = sqrt(diffVec(:,1).^2+diffVec(:,2).^2);
                    lengthFactor = repmat(lengthDiffVec,1,2);
                    stormPoints(idx2(firstLoc2),1:2)=stormPoints(idx2(firstLoc2),1:2)-(affectingFactor.*psfwidth-lengthFactor)/(affectingFactor.*psfwidth).*0.5.*diffVec;
                    stormPoints(idx2(secondLoc2),1:2)=stormPoints(idx2(secondLoc2),1:2)+(affectingFactor.*psfwidth-lengthFactor)/(affectingFactor.*psfwidth).*0.5.*diffVec;

                end

               stormPoints(idx2(firstLoc),:)=[];
                stormPoints = [stormPoints;meanCoords];
            end
        end
    end
    
end