function lengths = getFilamentLengths(lineCoords)
    lengths = [];
    for i = 1:size(lineCoords,2)
        lengths = [lengths,norm(lineCoords{i}(1,:) - lineCoords{i}(2,:))];
    end
end