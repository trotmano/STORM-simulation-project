function triangles = rescaleTriangles(triangles,scale)
    for i = 1:size(triangles,1)
        triangles{i} = triangles{i}*scale;
    end

end