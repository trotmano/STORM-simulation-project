function printAntibodies(listStartPoints, listEndPoints,fig, cAB)
    figure(fig)
    hold on
    for i = 1 : size(listStartPoints,1)
        line([ listStartPoints(i,1), listEndPoints(i,1)],...
                    [ listStartPoints(i,2),listEndPoints(i,2)],...
                    [ listStartPoints(i,3), listEndPoints(i,3)],'color',cAB)
    end
end