function objects = getCrossingLines(objects)
    for i = 1:2
        objects{size(objects,2)+1} = createLine();
    end
end

function line = createLine()
    for i = 1:2
        x = rand(1,1)*500;
        y = rand(1,1)*500;
        z = rand(1,1)*400;
        line(i,:) = [x,y,z];
    end
end