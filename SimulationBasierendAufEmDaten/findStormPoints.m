function [stormPoints, idxF, idxST] = findStormPoints(listEndPoints, abpf, sxy, sz,background)
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
        ilpmm3 = 30; %incorrect localizations per micrometer ^3
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
        intensities = abs(randn(size(listEndPoints,1),1))*3.8e4;
        listEndPoints = [listEndPoints;[x,y,z, intensities]];
    end
    
    nbrBlinkingEvents = randn(size(listEndPoints,1),1)*sqrt(abpf) + abpf;
    nbrBlinkingEvents(nbrBlinkingEvents<0) = 0;
    for i = 1:floor(max(max(nbrBlinkingEvents)))
        idx = nbrBlinkingEvents >= i;
        x = listEndPoints(idx,1) + randn(size(listEndPoints(idx),1),1) * sxy;
        y = listEndPoints(idx,2) + randn(size(listEndPoints(idx),1),1) * sxy;
        z = listEndPoints(idx,3) + randn(size(listEndPoints(idx),1),1) * sz;
        intensities = abs(randn(size(listEndPoints(idx),1),1))*3.8e4;
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
    
end