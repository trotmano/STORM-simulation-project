function writeStormPointsForVisp(stormPoints,fname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm)
    stormPoints = [stormPoints,zeros(size(stormPoints,1),1)];
    parts = strsplit(fname,'.');
    params = sprintf('loa-%3.3faoa-%3.3fbspnm-%3.3fpabs-%3.3fabpf-%3.3frof-%3.3fsxy-%3.3fsz-%3.3fbspsnm%3.3f',loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm);
    dlmwrite([parts{1},params,'_VispOutput.txt'],stormPoints,' ')
end