function printTriangles(triangles,fig,col)
    figure(fig)    
    hold on
    triangles = getMatrix(triangles);
    for i = 1:size(triangles,1)
        plot3(triangles(i,:,1),triangles(i,:,2),triangles(i,:,3),'color',col)
%                 hold on
    end
end