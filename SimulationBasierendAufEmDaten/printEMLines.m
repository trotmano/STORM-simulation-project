function printEMLines(points,fig,col)
    nbrPoints = size(points,2);
    figure(fig)
    hold on
    for i = 1:nbrPoints
        line(points{1,i}(:,1),points{1,i}(:,2),points{1,i}(:,3),'linewidth',2,'color',col)
    end
end