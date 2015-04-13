function triang = getMatrix(triangles)
%     triang = zeros(size(triangles,1),4,3);
    triang = zeros(size(triangles,1),3,3);
    for i = 1:size(triangles,1)
        tmp = triangles{i};
       triang(i,:,:) = tmp(1:3,1:3); 
    end
end