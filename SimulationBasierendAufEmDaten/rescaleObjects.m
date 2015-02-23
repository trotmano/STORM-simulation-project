function obj = rescaleObjects(obj,scale)
    for i = 1:size(obj,2)
        obj{i} = obj{i}*scale;
    end
end