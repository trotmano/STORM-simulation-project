function [basepoints,ep,idx] = findAntibodiesTri(triangles, bspsnm, pabs, loa, aoa, nocpsnm, docpsnm, areaClust, enableClust)

    triang = getMatrix(triangles);
%     triang = triang(1:1000,:,:);
    triang = triang(:,1:3,:);
    areas = getAreas(triang);
    if enableClust == 1
        [basepoints,idx,indicesNbrFluor,indicesPoints] = findBasePointsCluster(triang,areas,nocpsnm,docpsnm,areaClust);
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

function [points,indicesFluor,indicesNbrFluor,indicesPoints] = findBasePointsCluster(triang, areas, nocpsnm, docpsnm, areaClust)

    [basepointCluster,idIndices,idMatrix] = findRandomTriangles(triang,areas,nocpsnm);
    [neighMatrix,neighArea] = findNeighbouringTriangles(triang,idIndices,areaClust);
    [points,indicesFluor,indicesNbrFluor,indicesPoints] = setAntibodies(triang,idIndices,neighMatrix,neighArea,areas,docpsnm);
end

function [neighMatrix,neighArea] = findNeighbouringTriangles(triang,indicesClust_base,areaClust)

    lengthTriang = size(triang,1);
    lengthBase = size(indicesClust_base,2);
    neighMatrix = cell(lengthTriang,1);
    coordClust_base = zeros(lengthBase,3,3);
    for i = 1 : lengthBase
        coordClust_base(i,:,:) = triang(indicesClust_base(i),:,:);
    end
    
    triang_point1 = squeeze(triang(:,1,:));
    triang_point2 = squeeze(triang(:,2,:));
    triang_point3 = squeeze(triang(:,3,:));
    triang_allPoints = [triang_point1;triang_point2;triang_point3];
    for p = 1 : lengthBase
        coordClust_point = squeeze(coordClust_base(p,:,:));
        
        identicalPoints_point = ismember(triang_allPoints,coordClust_point,'rows');
        [row, ~ ] = find(identicalPoints_point);
        row = mod(row,lengthTriang);    
        
        identicalPoints = unique(row);
        triangId = indicesClust_base(p);
        [ ~ , deleteBase] = ismember(triangId,identicalPoints);
        identicalPoints(deleteBase) = 0;
        identicalPoints = nonzeros(identicalPoints);
        
        lenDiff = lengthTriang - length(identicalPoints);
        neighMatrix{triangId} = [identicalPoints; zeros(lenDiff,1)];
    end
    [neighMatrix,neighArea] = helper_RecursiveNeigh(triang,indicesClust_base,neighMatrix,areaClust);
end

function [neighMatrix,neighArea] = helper_RecursiveNeigh(triang,indicesClust_base, ...
                                                neighMatrix,areaClust)
    area = getAreas(triang);
    lengthBase = size(indicesClust_base,2);
    lengthTriang = size(triang,1);
    lengthMatrix = size(neighMatrix{indicesClust_base(1),1},1);
    neighArea = zeros(lengthTriang,1);
    for i = 1 : lengthBase
        i
        actual_area = 0;
        totalCount = 2;
        actual_neighborhood = nonzeros(neighMatrix{indicesClust_base(i),1}(:,totalCount - 1))';
        lenActualNeigh = length(actual_neighborhood);
        for j = 1 : lenActualNeigh
            actual_neighbor = actual_neighborhood(j);
            actual_area = actual_area + area(actual_neighbor);
        end
        
        while actual_area <= areaClust
            actual_neighborhood = nonzeros(neighMatrix{indicesClust_base(i),1}(:,totalCount - 1))';
            lenActualNeigh = length(actual_neighborhood);
            
            second_neighMatrix = helper_neighborhood(triang,actual_neighborhood);
            next_neighborhood = zeros(1,lengthTriang);
            begin = 1;
            for j = 1 : lenActualNeigh
                next_neighbors = nonzeros(second_neighMatrix{actual_neighborhood(j),1})';          
                diff_neighbors = setdiff(next_neighbors,next_neighborhood);
                lenDiff = length(diff_neighbors);
                eND = begin + lenDiff - 1;
                next_neighborhood(1,begin:eND) = diff_neighbors;
                begin = begin + lenDiff;
            end
            total_neighborhood = neighMatrix{indicesClust_base(i),1};
            diff_nextNeigh = setdiff(next_neighborhood,total_neighborhood)';
            lenDiff = lengthMatrix - length(diff_nextNeigh);
            newDiff_nextNeigh = [diff_nextNeigh ; zeros(lenDiff,1)];
            neighMatrix{indicesClust_base(i),1}(:,totalCount) = newDiff_nextNeigh;
            totalCount = totalCount + 1;
            for j = 1 : length(diff_nextNeigh)
                newActual_neighbor = diff_nextNeigh(j);
                actual_area = actual_area + area(newActual_neighbor);
            end
        end
        neighArea(indicesClust_base(i)) = actual_area;
    end
end

function neighMatrix = helper_neighborhood(triang,indicesClust_base)

    lengthTriang = size(triang,1);
    lengthBase = size(indicesClust_base,2);
    neighMatrix = cell(lengthTriang,1);
    coordClust_base = zeros(lengthBase,3,3);
    for i = 1 : lengthBase
        coordClust_base(i,:,:) = triang(indicesClust_base(i),:,:);
    end

    triang_point1 = squeeze(triang(:,1,:));
    triang_point2 = squeeze(triang(:,2,:));
    triang_point3 = squeeze(triang(:,3,:));
    triang_allPoints = [triang_point1;triang_point2;triang_point3];
    for p = 1 : lengthBase
        coordClust_point = squeeze(coordClust_base(p,:,:));

        identicalPoints_point = ismember(triang_allPoints,coordClust_point,'rows');
        [row, ~ ] = find(identicalPoints_point);
        row = mod(row,lengthTriang);

        identicalPoints = unique(row);
        triangId = indicesClust_base(p);
        [ ~ , deleteBase] = ismember(triangId,identicalPoints);
        identicalPoints(deleteBase) = 0;
        identicalPoints = nonzeros(identicalPoints);

        diff_length = lengthTriang - length(identicalPoints);
        neighMatrix{triangId} = [identicalPoints; zeros(diff_length,1)];
    end

%     lengthTriang = size(triang,1);
%     lengthBase = size(indicesClust_base,2);
%     neighMatrix = cell(lengthTriang,1);
%     coordClust_base = zeros(lengthBase,3,3);
%     for i = 1 : lengthBase
%         coordClust_base(i,:,:) = triang(indicesClust_base(i),:,:);
%     end
%     
%     totalCountNeighbors = 0;
%     for p = 1 : lengthBase
%         identical_coord = ismember(triang,coordClust_base(p,:,:));
%         lenCoord = length(identical_coord);
%         identitcalPoints = zeros(lenCoord,3);
%         for i = 1 : lenCoord
%             for j = 1 : 3
%                 identitcalPoints(i,j) = identical_coord(i,j,1)*identical_coord(i,j,2)*identical_coord(i,j,3);
%             end
%         end
%         sharedPoints = sum(identitcalPoints');
%         triangId = indicesClust_base(p);
%         neighMatrix{triangId} = zeros(lengthBase,1);
%         count = 1;
%         for k = 1 : length(sharedPoints)
%             if sharedPoints(k) > 0 && triangId ~= k
%                 neighMatrix{triangId}(count,1) = k;
%                 totalCountNeighbors = totalCountNeighbors +1;
%                 count = count + 1;
%             end
%         end
%     end
end

function [points,indices_Fluor,indices_nbrFluoph,indicesPoints] = setAntibodies(triang,indicesClust_base,neighMatrix,neighAllArea,areas,docpsnm)

    lengthBase = length(indicesClust_base);
    indices_nbrFluoph = zeros(length(triang),1);
    
    for i = 1 : lengthBase
        base = indicesClust_base(i);
        clusterArea = neighAllArea(base);
        neighborhood = nonzeros(neighMatrix{base,1});
        cluster_nbrFluoroph = floor(clusterArea * docpsnm);
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

    totalArea = sum(real(areas));
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

function neighMatrix = computeDistance(triang,indicesClust,neighMatrix)
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