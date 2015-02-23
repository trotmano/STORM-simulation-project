function points = importFilamentousStructures(fname)


fid = fopen(fname);

delimiter = ' ';
tline = fgets(fid);
points = {};
counter = 0;
while ischar(tline)
    parts = strsplit(tline,delimiter);
    if strcmp(parts(2),'#')&&strcmp(parts(3),'X')
        counter = counter + 1;
        tline = fgets(fid);
        parts = strsplit(tline,delimiter);
        tmp = [];
        while size(parts,2) == 6
            tmp = [tmp;[str2double(parts(3)),str2double(parts(4)),str2double(parts(5))]];
            tline = fgets(fid);
            parts = strsplit(tline,delimiter);
        end
        if (size(tmp)>0)
            points{counter,1} = tmp;
        else
            counter = counter -1;
        end
    end
    tline = fgets(fid);
end
points = points';
