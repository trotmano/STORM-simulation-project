function [basepoints,ep] = findFluorophoresTri(triangles, bspsnm, pabs, loa, aoa)
    ep = [];
    triang = getMatrix(triangles);
    areas = getAreas(triang);
    %areas = getAreas(triangles);
    nbrFluorophores = floor(sum(real(areas))*bspsnm*pabs);
    [basepoints,idx] = findBasePoints(nbrFluorophores,triang, areas);
    ep = getEndpoints(basepoints,triang,idx,loa, aoa);
end

function ep = getEndpoints(basepoints, triang,idx,loa, aoa)
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

function [points,idx] = findBasePoints(nbrFluorophores,triang,areas)
    
    points = zeros(nbrFluorophores,3);
    v1 = squeeze(triang(:,2,:) - triang(:,1,:));
    v2 = squeeze(triang(:,3,:) - triang(:,1,:));
    start = tic;
    %idx = getRandomTriangles2(areas, nbrFluorophores);
    toc(start)
    start = tic;
    idx = getRandomTriangles(areas, nbrFluorophores);
    toc(start)
%    idx = randi(size(triang,1),nbrFluorophores,1);
    for i = 1:nbrFluorophores
        while 1==1
            randx = rand(1);
            randy = rand(1);
            p = squeeze(triang(idx(i),1,:))+randx*v1(idx(i),:)'+randy*v2(idx(i),:)';
            A = [squeeze(triang(idx(i),1,:))';v1(idx(i),:);v2(idx(i),:)];
            x = linsolve(A',p);
            if x(2)<0||x(2)>1 ||x(3) <0||x(3)>1||x(2)+x(3)>1
%                 figure
%                 plot3(triang(idx(i),:,1),triang(idx(i),:,2),triang(idx(i),:,3))
%                 hold on
%                 plot3(p(1),p(2),p(3),'r*')
            else
                break
            end
        end
        points(i,:) = p';
    end
end

function idx = getRandomTriangles(areas, nbrFluorophores) %based on area of the triangles
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

function idx = getRandomTriangles2(areas, nbrFluorophores) %based on area of the triangles
    probs = areas/sum(areas);
    r = rand(nbrFluorophores,1);
    idx = arrayfun(@(z)sum(z >= cumsum([0, probs'])), r);
end


function areas = getAreas(triang)
    v1 = triang(:,2,:) - triang(:,1,:);
    v2 = triang(:,3,:) - triang(:,1,:);
    v2norm= normr(squeeze(v2));
    heigths = sqrt(sum(squeeze(v1).^2,2)-dot(squeeze(v1),v2norm,2).^2);
    areas = 0.5 * heigths .* sqrt(sum(squeeze(v2).^2,2));
end

function areas =  getAreas2(triangles)
    for i = 1:size(triangles)
        i
        currTriangle = triangles{i,1};
        v1 = currTriangle(2,:)-currTriangle(1,:);
        v2 = currTriangle(3,:)-currTriangle(1,:);
        height = sqrt(norm(v1)^2-(v1*(v2'/norm(v2)))^2);
        areas(i,1) = 0.5 * norm(v2)*height;
    end
end