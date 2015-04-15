function [basepoints,ep,idx] = findAntibodiesTri_3(triangles, bspsnm, pabs, loa, aoa, nocpsnm, doa, diameter, enableClust)

    triang = getMatrix(triangles);
%     triang = triang(1:5000,:,:);
    triang = triang(:,1:3,:);
    areas = getAreas(triang);
    if enableClust == 1
        [basepoints,idx,indicesNbrFluor,indicesPoints] = findBasePointsCluster(triang,areas,nocpsnm,doa,diameter);
        ep = getEndpointsClust(basepoints,triang,idx,indicesNbrFluor,indicesPoints,loa,aoa);
    else
        nbrFluorophores = floor(sum(real(areas))*bspsnm*pabs);
        [basepoints,idx] = findBasePoints(nbrFluorophores,triang,areas);
        ep = getEndpoints(basepoints,triang,idx,loa,aoa);
    end
end

function ep = getEndpoints(basepoints,triang,idx,loa, aoa)

    v1 = squeeze(triang(:,2,:) - triang(:,1,:));
    v2 = squeeze(triang(:,3,:) - triang(:,1,:));
    ep = zeros(size(basepoints,1),3);
    for i = 1:size(basepoints)
        vec = getVectorTri(aoa,loa);
        normTri = cross(v1(idx(i),:),v2(idx(i),:));
        finvec = findRotationTri(vec,normTri);
        ep(i,:) = basepoints(i,:)+finvec';
    end
end

function ep = getEndpointsClust(basepoints,triang,idx,idxNbr,idxPoints,loa,aoa)

    v1 = squeeze(triang(:,2,:) - triang(:,1,:));
    v2 = squeeze(triang(:,3,:) - triang(:,1,:));
    ep = zeros(size(basepoints,1),3);
    count = 1;
    for i = 1:size(idxNbr)
        nbrFluor = idxNbr(i);
        if nbrFluor == 0
            continue
        end
        for j = 1 : nbrFluor
            vec = getVectorTri(aoa,loa);
            normTri = cross(v1(i,:),v2(i,:));
            finvec = findRotationTri(vec,normTri);
            % weird counting procedure
            ep(count,:) = basepoints(count,:) + finvec';
            count = count + 1;
        end
    end
end

function [points,indicesFluor,indicesNbrFluor,indicesPoints] = findBasePointsCluster(triang, areas, nocpsnm, doa, diameter)

    [basepointCluster,indicesBaseClust,idMatrix] = findRandomTriangles(triang,areas,nocpsnm);
    clustEnvironment = buildClustEnvironment(triang,indicesBaseClust,diameter);
    [neighMatrix,neighArea] = findNeighbouringTriangles(triang,indicesBaseClust,clustEnvironment);
    [points,indicesFluor,indicesNbrFluor,indicesPoints] = setAntibodies(triang,indicesBaseClust,neighMatrix,neighArea,areas,doa);
end

function clustEnvironment = buildClustEnvironment(triang,indicesBaseClust,diameter)
    
    % Choose some random point within every triangle to
    % which the distance to base triangles of the clusters 
    % will be calculated
    pointInTriangles = squeeze(sum(triang(:,:,:),2))./3;
%     distance = sqrt(areaClust/pi);
    radius = diameter/2;
    clustEnvironment = computeDistance(triang,indicesBaseClust,pointInTriangles,radius);
end

function envMatrix = computeDistance(triang,indicesBaseClust,pointInTriangles,distance)
    % compute average distance of neighboring triangles
    % to base triangle
%     indicesClustBase = sort(indicesClust);
    lenBase = size(indicesBaseClust,2);
    lenTriang = size(triang,1);
    envMatrix = cell(lenBase,1);
    for i = 1 : lenBase
        baseTriangle = indicesBaseClust(i);
        baseCoord = pointInTriangles(baseTriangle,:);
        baseCoordMatrix = repmat(baseCoord,lenTriang,1);
        % fill environment Matrix with euclidean distance of every triangle
        % to the base triangle
        distanceMatrix = pointInTriangles - baseCoordMatrix;
        distanceMatrix = distanceMatrix.^2;
        distanceMatrix = sqrt(sum(distanceMatrix,2));
        tested_envTriang = find(distanceMatrix <= distance);
        envMatrix{i} = tested_envTriang;
    end
end

function [neighMatrix,neighArea] = findNeighbouringTriangles(triang,idBase_clust,clustEnvironment)

    lenTriang = size(triang,1);
    lenBase = size(idBase_clust,2);
    neighMatrix = cell(lenTriang,1);
    triangP1 = squeeze(triang(:,1,:));
    triangP2 = squeeze(triang(:,2,:));
    triangP3 = squeeze(triang(:,3,:));

    for p = 1 : lenBase
        % The coordinates of each base triangles will be compared to
        % the circumjacent triangles (clustEnvironment) which fulfill
        % the euclidean distance condition
        idBase = idBase_clust(p);
        coordBase = squeeze(triang(idBase,:,:));
        idEnv_base = clustEnvironment{p};
        lenEnv_base = size(idEnv_base,1);
        coordEnv_base = zeros(3*lenEnv_base,3);
        for i = 1 : lenEnv_base
            currentPoint = idEnv_base(i);
            coordEnv_base(i,:) = squeeze(triangP1(currentPoint,:,:));
            coordEnv_base(i + lenEnv_base,:) = squeeze(triangP2(currentPoint,:,:));
            coordEnv_base(i + 2*lenEnv_base,:) = squeeze(triangP3(currentPoint,:,:));
        end
        identPoints_logical = ismember(coordEnv_base,coordBase,'rows');
        [row, ~ ] = find(identPoints_logical);
        if size(row,1) == 0
            size(row,1)
            continue
        end
        row = mod(row,lenEnv_base);
        identicalPoints_1 = unique(row);
        if identicalPoints_1(1,1) == 0
            identicalPoints_1(1,1) = lenEnv_base;
        end
        
        lenIdentPoints = size(identicalPoints_1,1);
        identicalPoints = zeros(lenIdentPoints,1);
        for i = 1 : lenIdentPoints
            % Indices in environment to current base triangle will be
            % rewritten to actual indices of triangles
            current_idNeighbor = idEnv_base(identicalPoints_1(i));
            identicalPoints(i,1) = current_idNeighbor;
        end
        % delete base Triangle itself from neighborhood list
        [ ~ , deleteBase] = ismember(idBase,identicalPoints);
        if deleteBase ~= 0
            identicalPoints(deleteBase) = 0;
        else
            'Warning, base not included in env.'
        end
        
        identicalPoints = nonzeros(identicalPoints);
        % write neighborhood list in a cell (neighMatrix), which is accesible
        % by the indice of the base triangle of a cluster
        
        %         lenDiff = lenTriang - length(identicalPoints);
        lenClustEnv = length(idEnv_base);
        extIdPoints = zeros(lenClustEnv,1);
        for i = 1 : length(identicalPoints)
            extIdPoints(i) = identicalPoints(i);
        end
        neighMatrix{idBase} = extIdPoints;
        %         neighMatrix{idBase} = [identicalPoints; zeros(lenDiff,1)];
    end
    [neighMatrix,neighArea] = helper_RecursiveNeigh(triang,idBase_clust,clustEnvironment, ...
                                                    neighMatrix,triangP1,triangP2,triangP3);
end

function [neighMatrix,neighArea] = helper_RecursiveNeigh(triang,idBase_clust,clustEnvironment, ...
                                                        neighMatrix,triP1,triP2,triP3)
    area = getAreas(triang);
    lenBase = size(idBase_clust,2);
    lenTriang = size(triang,1);
    neighArea = zeros(lenTriang,1);
    for i = 1 : lenBase
        i
        newRow = 2;
        actual_area = 0;
        current_neighborhood = nonzeros(neighMatrix{idBase_clust(i),1}(:,newRow - 1))';
        lenActualNeigh = length(current_neighborhood);
        for j = 1 : lenActualNeigh
            actual_neighbor = current_neighborhood(j);
            actual_area = actual_area + area(actual_neighbor);
        end
        
        while true
            current_neighborhood = nonzeros(neighMatrix{idBase_clust(i),1}(:,newRow - 1))';
            lenCurrentNeigh = length(current_neighborhood);
            next_neighborhood = zeros(1,lenTriang);
            
            currentBase = i;
            second_neighMatrix = helper_neighborhood2(triang,current_neighborhood,clustEnvironment, ...
                                                        currentBase,triP1,triP2,triP3);
            begin = 1;
            for j = 1 : lenCurrentNeigh
                next_neighbors = nonzeros(second_neighMatrix{current_neighborhood(j),1})';          
                diff_neighbors = setdiff(next_neighbors,next_neighborhood);
                lenDiff = length(diff_neighbors);
                eND = begin + lenDiff - 1;
                next_neighborhood(1,begin:eND) = diff_neighbors;
                begin = begin + lenDiff;
            end
            total_neighborhood = neighMatrix{idBase_clust(i),1};
            diff_nextNeigh = setdiff(next_neighborhood,total_neighborhood)';
            if size(diff_nextNeigh,1) == 0
                break
            end
            
            %             lenDiff = lenMatrix - length(diff_nextNeigh);
            %             newDiff_nextNeigh = [diff_nextNeigh ; zeros(lenDiff,1)];
            %             neighMatrix{idBase_clust(i),1}(:,newRow) = newDiff_nextNeigh;
            
            lenClustEnv = size(neighMatrix{idBase_clust(i),1},1);
            ext_diffNextNeigh = zeros(lenClustEnv,1);
            for j = 1 : length(diff_nextNeigh)
                ext_diffNextNeigh(j) = diff_nextNeigh(j);
            end
            neighMatrix{idBase_clust(i),1}(:,newRow) = ext_diffNextNeigh;
            newRow = newRow + 1;
            
            for j = 1 : length(diff_nextNeigh)
                newActual_neighbor = diff_nextNeigh(j);
                actual_area = actual_area + area(newActual_neighbor);
            end
        end
        neighArea(idBase_clust(i)) = actual_area;
    end
end

function neighMatrix = helper_neighborhood(triang,idBase_clust,clustEnvironment,currentBase, ...
                                             triP1,triP2,triP3)

    lenTriang = size(triang,1);
    lenBase = size(idBase_clust,2);
    neighMatrix = cell(lenTriang,1);
    idEnv_base = clustEnvironment{currentBase};

    for p = 1 : lenBase
        % The coordinates of each base triangles will be compared to
        % the circumjacent triangles (clustEnvironment) which fulfill
        % the euclidean distance condition
        idBase = idBase_clust(p);
        coordBase = squeeze(triang(idBase,:,:));
        lenEnv_base = size(idEnv_base,1);
        coordEnv_base = zeros(3*lenEnv_base,3);
        for i = 1 : lenEnv_base
            currentPoint = idEnv_base(i);
            coordEnv_base(i,:) = squeeze(triP1(currentPoint,:,:));
            coordEnv_base(i + lenEnv_base,:) = squeeze(triP2(currentPoint,:,:));
            coordEnv_base(i + 2*lenEnv_base,:) = squeeze(triP3(currentPoint,:,:));
        end
        identPoints_logical = ismember(coordEnv_base,coordBase,'rows');
        [row, ~ ] = find(identPoints_logical);
        if size(row,1) == 0
            size(row,1)
            continue
        end
        row = mod(row,lenEnv_base);
        identicalPoints_1 = unique(row);
        if identicalPoints_1(1,1) == 0
            identicalPoints_1(1,1) = lenEnv_base;
        end
        
        lenIdentPoints = size(identicalPoints_1,1);
        identicalPoints = zeros(lenIdentPoints,1);
        for i = 1 : lenIdentPoints
            % Indices in environment to current base triangle will be
            % rewritten to actual indices of triangles
            current_idNeighbor = idEnv_base(identicalPoints_1(i));
            identicalPoints(i,1) = current_idNeighbor;
        end
        % delete base Triangle itself from neighborhood list
        [ ~ , deleteBase] = ismember(idBase,identicalPoints);
        if deleteBase ~= 0
            identicalPoints(deleteBase) = 0;
        else
            'Warning, base not included in env.'
        end
        
        identicalPoints = nonzeros(identicalPoints);
        % write neighborhood list in a cell (neighMatrix), which is accesible
        % by the indice of the base triangle of a cluster
        
        lenDiff = lenTriang - length(identicalPoints);
        neighMatrix{idBase} = [identicalPoints; zeros(lenDiff,1)];
    end
end

function neighMatrix = helper_neighborhood2(triang,idBase_clust,clustEnvironment,currentBase, ...
                                                triP1,triP2,triP3)

    lenTriang = size(triang,1);
    lenBase = size(idBase_clust,2);
    neighMatrix = cell(lenTriang,1);
    % Generate list containing coordinates for all
    % elements in cluster environment
    idEnv_base = clustEnvironment{currentBase};
    lenEnv_base = size(idEnv_base,1);
    coordEnv_base = zeros(3*lenEnv_base,3);
    for i = 1 : lenEnv_base
        currentPoint = idEnv_base(i);
        coordEnv_base(i,:) = squeeze(triP1(currentPoint,:,:));
        coordEnv_base(i + lenEnv_base,:) = squeeze(triP2(currentPoint,:,:));
        coordEnv_base(i + 2*lenEnv_base,:) = squeeze(triP3(currentPoint,:,:));
    end

    for p = 1 : lenBase
        % The coordinates of each base triangles will be compared to
        % the circumjacent triangles (clustEnvironment) which fulfill
        % the euclidean distance condition
        idBase = idBase_clust(p);
        coordBase = squeeze(triang(idBase,:,:));
        identPoints_logical = ismember(coordEnv_base,coordBase,'rows');
        [row, ~ ] = find(identPoints_logical);
        if size(row,1) == 0
            size(row,1)
            continue
        end
        row = mod(row,lenEnv_base);
        identicalPoints_1 = unique(row);
        if identicalPoints_1(1,1) == 0
            identicalPoints_1(1,1) = lenEnv_base;
        end
        
        lenIdentPoints = size(identicalPoints_1,1);
        identicalPoints = zeros(lenIdentPoints,1);
        for i = 1 : lenIdentPoints
            % Indices in environment to current base triangle will be
            % rewritten to actual indices of triangles
            current_idNeighbor = idEnv_base(identicalPoints_1(i));
            identicalPoints(i,1) = current_idNeighbor;
        end
        % delete base Triangle itself from neighborhood list
        [ ~ , deleteBase] = ismember(idBase,identicalPoints);
        if deleteBase ~= 0
            identicalPoints(deleteBase) = 0;
        else
            'Warning, base not included in env.'
        end
        
        identicalPoints = nonzeros(identicalPoints);
        % write neighborhood list in a cell (neighMatrix), which is accesible
        % by the indice of the base triangle of a cluster
        
        %         lenDiff = lenTriang - length(identicalPoints);
        lenClustEnv = length(idEnv_base);
        extIdPoints = zeros(lenClustEnv,1);
        for i = 1 : length(identicalPoints)
            extIdPoints(i) = identicalPoints(i);
        end
        neighMatrix{idBase} = extIdPoints;
        %         neighMatrix{idBase} = [identicalPoints; zeros(lenDiff,1)];                                                                                                                                                                                                                                                                                                                                     
    end
end

% function envMatrix = computeAverageDistance(triang,indicesBaseClust,pointInTriangles,distance)
%     % compute average distance of neighboring triangles
%     % to base triangle
% %     indicesClustBase = sort(indicesClust);
%     lenBase = size(indicesBaseClust,2);
%     lenTriang = size(triang,1);
%     envMatrix = cell(lenBase,1);
%     for i = 1 : lenBase
%         baseTriangle = indicesBaseClust(i);
%         baseCoord = pointInTriangles(baseTriangle,:);
%         baseCoordMatrix = repmat(baseCoord,lenTriang,1);
%         % fill environment Matrix with euclidean distance of every triangle
%         % to the base triangle
%         distanceMatrix = pointInTriangles - baseCoordMatrix;
%         distanceMatrix = distanceMatrix.^2;
%         distanceMatrix = sqrt(sum(distanceMatrix,2));
%         tested_envTriang = find(distanceMatrix <= distance);
%         envMatrix{i} = tested_envTriang;
%     end
% end

function [points,indices_Fluor,indices_nbrFluoph,indicesPoints] = setAntibodies(triang,indicesClust_base,neighMatrix,neighAllArea,areas,doa)

    lengthBase = length(indicesClust_base);
    indices_nbrFluoph = zeros(length(triang),1);
    
    for i = 1 : lengthBase
        base = indicesClust_base(i);
        clusterArea = neighAllArea(base);
        neighborhood = nonzeros(neighMatrix{base,1});
        cluster_nbrFluoroph = floor(clusterArea * doa);
        cluster_nbrFluoroph
        lenNeigh = length(neighborhood);
        
        startingindices = [];
        startingsum = [];
        parts = 10;
        partsum = 0;
        counter = 1;
        for j = 1 : lenNeigh
            actual_neigh = neighborhood(j);
            %             actual_neighArea = areas(actual_neigh);
            %             fraction = actual_neighArea / clusterArea;
            partsum = partsum + areas(actual_neigh);
            divider = clusterArea / parts * counter;
            if partsum > divider
                startingindices(counter) = j;
                startingsum(counter) = partsum;
                counter = counter + 1;
            end
        end
        randd = rand(cluster_nbrFluoroph,1);
        randd = randd .* clusterArea;
        sizeStartingSum = size(startingsum,2);
        for j = 1 : cluster_nbrFluoroph
            randD = randd(j);
            startId = 1;
            for k = 2 : sizeStartingSum
               if randD > startingsum(k)
                   startId = k-1;
               end
            end
            partsum = startingsum(startId);
            for k = startingindices(startId) : lenNeigh
                neighbor = neighborhood(k);
                partsum = partsum + areas(neighbor);
                if randD < partsum
                    % Intersection of clusters raises fluorophore number
                    indices_nbrFluoph(neighbor) = indices_nbrFluoph(neighbor) + 1;
                    break
                end
            end
        end
    end
    
    nbrFluorophores = sum(indices_nbrFluoph);
    v1 = squeeze(triang(:,2,:) - triang(:,1,:));
    v2 = squeeze(triang(:,3,:) - triang(:,1,:));
    points = zeros(nbrFluorophores,3);
    indicesPoints = zeros(nbrFluorophores,1);
    count = 1;
    for i = 1 : length(indices_nbrFluoph)
        actual_nbrFluor = indices_nbrFluoph(i);
        if actual_nbrFluor == 0
            continue
        end
        for j = 1 : actual_nbrFluor
            randx = rand(1);
            randy = rand(1);
            while randx + randy >= 1
                randx = rand(1);
                randy = rand(1);
            end
            p = squeeze(triang(i,1,:))+randx*v1(i,:)'+randy*v2(i,:)';
            points(count,:) = p';
            indicesPoints(count) = i;
            count = count + 1;
        end
    end
    [row, ~ ] = find(indices_nbrFluoph);
    indices_Fluor = row;
end

function [points,idIndices,idMatrix] = findRandomTriangles(triang,areas,nocpsnm)
%     rng('shuffle');
    totalArea = sum(real(areas))
    numClusters = floor(nocpsnm * totalArea);   %number of clusters;
    len = length(triang(:,1,1));
    while numClusters == 0
        'number of Clusters is zero, as nocpsnm is too small!'
        nocpsnmInput = input('Type in a higher nocpsnm parameter: ');
        numClusters = floor(nocpsnmInput * totalArea);
    end
    if numClusters > len
        'Warning: All Triangles are part of a cluster!'
        numClusters = len;
    end
    'number of Clusters: '
    numClusters
    idIndices = randperm(len,numClusters);
    idIndices = sort(idIndices);
    idMatrix = zeros(numClusters,3,3);
    for i = 1 : numClusters
        k = idIndices(i);
        idMatrix(i,:,:) = triang(k,:,:);
    end
    v1 = squeeze(triang(:,2,:) - triang(:,1,:));
    v2 = squeeze(triang(:,3,:) - triang(:,1,:));
    points = zeros(numClusters,3);
    for i = 1:numClusters
        randx = rand(1);
        randy = rand(1);
        while randx + randy >= 1
            randx = rand(1);
            randy = rand(1);
        end
        p = squeeze(triang(idIndices(i),1,:))+randx*v1(idIndices(i),:)'+randy*v2(idIndices(i),:)';
        points(i,:) = p';
    end
end

function idx = getRandomTriangles(areas, nbrFluorophores)
    tot = sum(areas);
    idx = [];
    startingindices = [];
    startingsum = [];
    partsum = 0;
    counter = 1;
    parts = 10000;
    for i = 1:size(areas)
        partsum = partsum + areas(i);
        if partsum > tot/parts*counter
            startingindices(counter) = i;
            startingsum(counter) = partsum;
            counter = counter + 1;
        end
    end
    randd = rand(nbrFluorophores,1);
    for i = 1:nbrFluorophores
        randD = randd(i)*tot;
        for k = 2:size(startingsum,2)
            if randD>startingsum(k)
                startidx = k-1;
            end
        end
        partsum = startingsum(startidx);
        for j = startingindices(startidx):size(areas)
            partsum = partsum + areas(j);
            if randD<partsum
                idx = [idx ,j-1];
                break
            end
        end
    end
end

function [points,idx] = findBasePoints(nbrFluorophores,triang,areas)
    v1 = squeeze(triang(:,2,:) - triang(:,1,:));
    v2 = squeeze(triang(:,3,:) - triang(:,1,:));
    start = tic;
    %idx = getRandomTriangles2(areas, nbrFluorophores);
    toc(start)
    start = tic;
    idx = getRandomTriangles(areas, nbrFluorophores);
    toc(start)
%    idx = randi(size(triang,1),nbrFluorophores,1);
    points = zeros(nbrFluorophores,3);
    for i = 1:nbrFluorophores
        randx = rand(1);
        randy = rand(1);
        while randx + randy >= 1
            randx = rand(1);
            randy = rand(1);
        end
        p = squeeze(triang(idx(i),1,:))+randx*v1(idx(i),:)'+randy*v2(idx(i),:)';
%         A = [squeeze(triang(idx(i),1,:))';v1(idx(i),:);v2(idx(i),:)];
%         x = linsolve(A',p);
        points(i,:) = p';
    end
end

function areas = getAreas(triang)
    v1 = triang(:,2,:) - triang(:,1,:);
    v2 = triang(:,3,:) - triang(:,1,:);
    v2norm = normr(squeeze(v2));
    heigths = sqrt(sum(squeeze(v1).^2,2)-dot(squeeze(v1),v2norm,2).^2);
    areas = 0.5 * heigths .* sqrt(sum(squeeze(v2).^2,2));
end

function neighMatrix = computeDistance2(triang,indicesClust,neighMatrix)
    % compute average distance of neighboring triangles
    % to base triangle
    indicesClustBase = sort(indicesClust);
    lengthBase = size(indicesClust,2);
    lengthNeighMatrix = size(neighMatrix{indicesClustBase(1),1},1);
    for i = 1 : lengthBase
        count = 0;
        baseTriangle = indicesClustBase(i);
        baseCoord = neighMatrix{baseTriangle};
        for j = 1 : lengthNeighMatrix
            neighTriangle = baseCoord(j,1);
            if neighTriangle == 0
                break
            end
            coordBaseTriang = squeeze(mean(triang(baseTriangle,:,:)));
            coordNeighTriang = squeeze(mean(triang(neighTriangle,:,:)));
            distance = pdist([coordBaseTriang,coordNeighTriang]','euclidean');
            neighMatrix{baseTriangle}(j,2) = distance;
            count = count + 1;
        end
        sumDistances = sum(neighMatrix{baseTriangle}(:,2));
        neighMatrix{baseTriangle}(1,3) = sumDistances/count;
    end
end