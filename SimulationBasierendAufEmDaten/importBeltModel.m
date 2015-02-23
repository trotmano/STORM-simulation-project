function points2 = importBeltModel(fname)


fid = fopen(fname);

delimiter = ' ';
tline = fgets(fid);
counter = 0;
points = {};
tmp = [];
tmpPunktScatter = [];
while ischar(tline)
    counter = counter +1;
    tmp2 = strsplit(tline,' ');
    tmpPunktScatter(1,1) = str2double(tmp2{1,3});
    tmpPunktScatter(1,2) = str2double(tmp2{1,4});
    tmpPunktScatter(1,3) = str2double(tmp2{1,5});
    tline = fgets(fid);
    if tline == -1
        break
    end
    tmp2 = strsplit(tline,delimiter);
    tmpPunktScatter(2,1) = str2double(tmp2{1,3});
    tmpPunktScatter(2,2) = str2double(tmp2{1,4});
    tmpPunktScatter(2,3) = str2double(tmp2{1,5});
    if tmpPunktScatter(1,3)==tmpPunktScatter(2,3)
        points{size(points,2)+1} = tmpPunktScatter;
    else
        
    end
end
points2 = {};
for i = 1:size(points,2)
    points2{size(points2,2)+1} = points{i};
    points2{size(points2,2)+1} = points{i}+[0,0,-4;0,0,-1];
end
fclose(fid);

