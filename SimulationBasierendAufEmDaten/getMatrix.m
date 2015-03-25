function triang = getMatrix(triangles)
    triang = zeros(size(triangles,1),4,3);
%     triang = zeros(size(triangles,1),3,3);
    for i = 1:size(triangles,1)
       triang(i,:,:) = triangles{i}; 
    end
end