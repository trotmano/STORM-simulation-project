function [length, cummulativeLengths] = getLengthOfStructure(lines)
    for i = 1:size(lines,1)-1
        if i > 1
            cummulativeLengths(i) = norm([lines(i,1) - lines(i+1,1),...
                                         lines(i,2) - lines(i+1,2),...
                                         lines(i,3) - lines(i+1,3)]) + cummulativeLengths(i-1);
        else
            cummulativeLengths(i) = norm([lines(i,1) - lines(i+1,1),...
                                         lines(i,2) - lines(i+1,2),...
                                         lines(i,3) - lines(i+1,3)]);
            length = cummulativeLengths(1);
        end
    end
    length = cummulativeLengths(end);
end