function triangles = swapColumns(triangles,c1,c2)
    for i = 1:size(triangles,1)
        tmp = triangles{i}(:,c1);
        triangles{i}(:,c1) = triangles{i}(:,c2);
        triangles{i}(:,c2) = tmp;
    end

end