function points = importTriangles(fname)
basename = strsplit(fname,'.');
ext = basename{2};
basename = basename{1};
if exist([basename,'_rs.mat'])
    points = load([basename,'_rs.mat']);
    points = points.points;
else
    fid = fopen(fname);

    delimiter = ' ';
    tline = fgets(fid);
    %tline = fgets(fid);
    %tline = fgets(fid);
    %tline = fgets(fid);
    counter = 0;
    points = {};
    tmp = [];
    tmpPunktScatter = [];
    tmplines = [];
    while ischar(tline)
        if findstr(tline,'pp 3')
            counter = counter + 1;
            tline = fgets(fid);
            for i = 1:3
                tmp = strsplit(tline,delimiter);
                tmplines = [tmplines;tmp(1:3)];
                tline = fgets(fid);
            end
            tmplines = [tmplines;tmplines(1,:)];
            tmplines = str2double(tmplines);
            if(sum(tmplines(2,:) == tmplines(3,:))==3||sum(tmplines(1,:) == tmplines(3,:))||sum(tmplines(2,:) == tmplines(1,:))) %all three coordinates are equal
                size(points)
                counter = counter -1;
            else
                points{counter,1} = tmplines;
            end
            
        else
            tline = fgets(fid);
        end
       tmplines = [];
       counter
    end
    save([basename,'_rs.mat'],'points')
end



end