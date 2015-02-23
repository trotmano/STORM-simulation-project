function writeOutputFileMalk(stormPoints,fname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm)
    stormPoints = [stormPoints(:,1:3),zeros(size(stormPoints,1),1),stormPoints(:,4)];
    parts = strsplit(fname,'.'); 
    params = sprintf('loa-%3.3faoa-%3.3fbspnm-%3.3fpabs-%3.3fabpf-%3.3frof-%3.3fsxy-%3.3fsz-%3.3fbspsnm%3.3f',loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm);
    outputFname = [parts{1},params,'_MalkOutput.txt'];
    fid = fopen(outputFname,'w+t');
    fprintf(fid,'# <localizations insequence="true" repetitions="variable"><field identifier="Position-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in X" unit="nanometer" min="0 m" max="2.2078e-005 m" /><field identifier="Position-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in Y" unit="nanometer" min="0 m" max="1.8487e-005 m" /><field identifier="Position-2-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in Z" unit="nanometer" min="0 m" max="9.1e-007 m" /><field identifier="ImageNumber-0-0" syntax="integer" semantic="frame number" unit="frame" min="0 fr" /><field identifier="Amplitude-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="emission strength" unit="A/D count" /></localizations>');
    fprintf(fid,'\n');
    fclose(fid);
    dlmwrite(outputFname,stormPoints,'-append','Delimiter',' '); 
end


