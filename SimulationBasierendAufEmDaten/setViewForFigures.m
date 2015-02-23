function setViewForFigures(az,el)
    figs = findobj('Type','figure');
    for i = 1:size(figs,2)
        figure(figs(i))
        view([az,el])
        axis equal
        %axis([600-40,800+40,650-40,1000+40,0-40,120+40])
    end
end