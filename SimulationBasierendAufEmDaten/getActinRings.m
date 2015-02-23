function objects = getActinRings(objects)
    for i = 1:10
        objects{size(objects,2)+1} = createRing(0+i*190,0,0,100,300);
    end
end

function ring = createRing(x,y,z,r,n)
    deltaPhi = 2*pi / n; %increment in angle 
    
    figure(1)
    hold on
    for i = 1:n+1 %n+1 to close the circle
        tmp = [x,y+r*cos(deltaPhi*i),z+r*sin(deltaPhi*i)];
        ring(i,:) = tmp;
    end

end