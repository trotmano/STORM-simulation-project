function obj = rescaleObjects(obj,scale)
    for i = 1:size(obj,1)
        obj{i} = obj{i}*scale;
    end
end