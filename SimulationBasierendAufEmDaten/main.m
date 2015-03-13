function main()
    
    
	%rng(5);
    loa = 1; %length of both antibodies combined
    aoa = 90/180*pi; %angle of antibody
    %bspnm = 1.65; %binding sites per nm  %for microtubules with only alpha tubuli stained it should be like 13/8 1.65
    bspnm = 1/2.75;
    pabs = 0.51; %part of available binding sites
    abpf = 14;% average blinking per fluorophor
    rof = 3.5;%radius of filament
    fpab = 1.5; %fluorophores per antibody

    cEM = [0,0,0]; %color of EM data
    cSTORM = [1,0,0]; %color of STORM result
    cAB = [0,1,0]; %color of Antibody
    sxy = 1; %sigma of fitting error in xy direction
    sz = 1; %sigma of fitting error in z direction
    
    doc = 0; %degree of clustering, part of all localizations that are clustered
    nocpsmm = 1; %number of clusters per square micrometer
    docpsnm = 0.01; %denstiy of clusters in antibodies per square nm
    bd = 50/10000/10000; %blinking density in number fluorophores per square nm
    
    bspsnm = 10/600.;%.0159/2; %binding sites per square nanometer
    fname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Synaptic Actin Data- Experimental + Simulation\EM Models\Newest Model- 150209\150209-nff';
    outputname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Actin in Infected Erythrocytes\ExtractedSurface-Different Outputs\output.txt';
    %objects = importTriangles(fname);
    %objects = getActinRings({});
    %objects = getCrossingLines({});
    objects = importEMData(fname);
    %objects = importPly(fname);
    %objects = importFilamentousStructures(fname);
    %objects = objects{1};
    %objects = swapColumns(objects,2,3);
    %objects = rescaleObjects(objects,2.3); %WICHTIG FÜR MITOCHONDRIEN
    %objects = importSelfMeassuredMicrotubuli();
   
    [ap,ep,stormPoints] = doSimulation(objects, bspsnm, pabs,loa,aoa,doc,nocpsmm,docpsnm,bspnm,rof,abpf,sxy,sz,bd,fpab);
    if (size(stormPoints,1)>0)
        writeOutput(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    end
    visualizeResults(objects,ap,ep,stormPoints)
end

function visualizeResults(objects,ap,ep,stormPoints)
    figHandle = findobj('Type','figure');
    close(figHandle);
    showEM(objects)
    showAntibodies(ap,ep)
    showEMAntibodies(objects,ap,ep)
    showStormPoints(stormPoints)
    showEMAntibodiesStormPoints(objects,ap,ep,stormPoints)

    setViewForFigures(150,30)

end

function writeOutput(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab)
    writeStormPointsFRC(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    writeStormPointsForVisp(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    writeOutputFileMalk(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    writeStormPointsForAmira(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
end

function [ap,ep,stormPoints] = doSimulation(objects, bspsnm, pabs,loa,aoa,doc,nocpsmm,docpsnm,bspnm,rof,abpf,sxy,sz,bd,fpab)
    if isSurfaceData(objects)
        [ap,ep] = findAntibodiesTri(objects, bspsnm, pabs, loa, aoa, doc, nocpsmm, docpsnm);
    else
        [ap,ep,idx]=findLines(objects, bspnm, pabs,aoa,loa,rof);
    end
    [stormPoints, idxF ,idxSt] = findStormPoints(ep, abpf, sxy, sz, bd,fpab,true);
end

function showEMAntibodiesStormPoints(objects,ap,ep,stormPoints)
    fig = figure(findNextFigure);
    clf
    hold on
    if isSurfaceData(objects)
        printTriangles(objects,fig)
    else
        printEMLines(objects,fig,[0,0,0])
    end
    printAntibodies(ap, ep,fig, [0,1,0])
    printSTORMPoints(stormPoints,fig)
end

function showStormPoints(stormPoints)
    fig = figure(findNextFigure);
    clf
    hold on
    printSTORMPoints(stormPoints,fig)
end

function showEMAntibodies(objects,ap,ep)
    fig = figure(findNextFigure);
    clf
    hold on
    if isSurfaceData(objects)
        printTriangles(objects,fig,[0,0,0])
    else
        printEMLines(objects,fig,[0,0,0])
    end
    printAntibodies(ap, ep,fig, [0,1,0])
end

function showAntibodies(ap,ep)
    fig = figure(findNextFigure);
    clf
    hold on
    printAntibodies(ap, ep,fig, [0,1,0])
end

function showEM(objects)
    fig = figure(findNextFigure);
    clf
    hold on
    if isSurfaceData(objects)
        printTriangles(objects,fig,[0,0,0])
    else
        printEMLines(objects,fig,[0,0,0])
    end
end

function issd = isSurfaceData(objects)
    issd = 1;
    for i = 1:size(objects,2)
        if (size(objects{i},1)>4 || size(objects{i},1)==2)
            issd = 0;
            break;
        end
    end
end

function idx = findNextFigure()
    figureHandles = findobj('Type','figure');
    if size(figureHandles,1) == 0
        idx = 1;
    else
        idx = max(figureHandles)+1;
    end
end
