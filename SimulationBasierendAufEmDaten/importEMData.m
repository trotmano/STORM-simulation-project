function points = importEMData(fname)

saveFile = 1;
if exist([fname,'_rs'], 'file') == 2
    fname = [fname,'_rs'];
    saveFile = 0;
end
fid = fopen(fname);

delimiter = ' ';
tline = fgets(fid);
points = {};
tmp = [];
tmpPunktScatter = [];
counter = 0;
while ischar(tline)
    if size(tline,2)>=3 && ~isempty(findstr(tline(1:3), 'p 3'))
        tmp = [];
        tline = fgets(fid);
        tmp(size(tmp,1)+1,:) = str2double(strsplit(tline,delimiter));
        tline = fgets(fid);
        tmp(size(tmp,1)+1,:) = str2double(strsplit(tline,delimiter));
        tline = fgets(fid);
        %tmp(size(tmp,1)+1,:) = str2double(strsplit(tline,delimiter));
        points{size(points,2)+1} = tmp;
    end
    
    if size(tline,2)>=2 && ~isempty(findstr(tline(1:2), 's '))
        tmp2 = strsplit(tline,delimiter);
        tmpPunktScatter(1,1) = str2double(tmp2{1,2});
        tmpPunktScatter(1,2) = str2double(tmp2{1,3});
        tmpPunktScatter(1,3) = str2double(tmp2{1,4});
        tline = fgets(fid);
        tmp2 = strsplit(tline,delimiter);
        tmpPunktScatter(2,1) = str2double(tmp2{1,2});
        tmpPunktScatter(2,2) = str2double(tmp2{1,3});
        tmpPunktScatter(2,3) = str2double(tmp2{1,4});
        points{size(points,2)+1} = tmpPunktScatter;
    end
   
    tline = fgets(fid);
end

fclose(fid);

if saveFile
    fid = fopen([fname,'_rs'],'w');
    for i = 1:size(points,2)
        if size(points{i},1) == 3
            fprintf(fid, 'p 3\n');
            for j = 1:3
                fprintf(fid, '%4.2f %4.2f %4.2f\n', points{1,i}(j,1), points{1,i}(j,2), points{1,i}(j,3));
            end
        else
            fprintf(fid, 'scatter\n');
            for j = 1:2
                fprintf(fid, 's %4.2f %4.2f %4.2f\n', points{1,i}(j,1), points{1,i}(j,2), points{1,i}(j,3));
            end
        end
    end
    fclose(fid);
end