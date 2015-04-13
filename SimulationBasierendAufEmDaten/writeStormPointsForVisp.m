function writeStormPointsForVisp(stormPoints,fname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab)
    %stormPoints = [stormPoints,zeros(size(stormPoints,1),1)];
    parts = strsplit(fname,'.');
    params = sprintf('loa%2.2faoa%2.2fbspnm%2.2fpabs%2.2fabpf%2.2frof%2.2fsxy%2.2fsz%2.2fbspsnm%2.2ffpab%2.2f',loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    dlmwrite([parts{1},params,'_VispOutput.txt'],stormPoints,' ')
end