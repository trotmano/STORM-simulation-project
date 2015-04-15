function main()
profile on
    figHandle = findobj('Type','figure');
    close(figHandle);
% 	rng(5);
    rng('shuffle');
    
    loa = 10; %length of both antibodies combined
    aoa = 90/180*pi; %angle of antibody
    bspnm = 1.65; %binding sites per square nanometer -- for microtubules with only alpha tubuli stained it should be like 13/8 = 1.65
    %bspnm = 1/2.75;
    %bspsnm = .0159/2;
    pabs = 0.51; %part of available binding sites
    abpf = 14;% average blinking per fluorophor
    rof = 3.5;%radius of filament
    fpab = 1.5; %fluorophores per antibody

    cEM = [0,0,0]; %color of EM data
    cSTORM = [1,0,0]; %color of STORM result
    cAB = [0,1,0]; %color of Antibody
    
    sxy = 8; %sigma of fitting error in xy direction
    sz = 35; %sigma of fitting error in z direction
    
    enableClust = 1;
    noc = 0.0001; %number of clusters per square nm
    %areaClust = 3900; %average area of a cluster
    diaClust = 70; %average diameter of a cluster
    doa = 0.05; %denstiy of antibodies in clusters per square nm
        
    bd = 50/10000/10000; %blinking density in number fluorophores per square nm
    bspsnm = 10/600.;%.0159/2; %binding sites per square nanometer
    
    %fname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Organelle Library\Microtubules\Microtubules.wimp';
%     fname = '/media/Dev_d/Persönlicher Ordner/Docs/Skripte/Master/Studium/S_01/Kuner/EM Tomography Model/Mitochondria-Tomogram-beta-islet-cells.nff';
    fname = '/media/Dev_d/Persönlicher Ordner/Docs/Skripte/Master/Studium/S_01/Kuner/EM Tomography Model/Mitochondria-Tomogram-beta-islet-cells-scaled10000.nff';
%     fname = '/media/Dev_d/Persönlicher Ordner/Docs/Skripte/Master/Studium/S_01/Kuner/EM Tomography Model/mitoBig.nff';

    %outputname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Organelle Library\Mitochondria\STORM Simulation\Mitochondria-Tomogram-beta-islet-cells.nff';
    outputname = '/media/Dev_d/Persönlicher Ordner/Docs/Skripte/Master/Studium/S_01/Kuner/EM Tomography Model/Output/mitoscaled1_noc0.0001_diaClust70_doa0.05.nff';
    
    %fname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Synaptic Actin Data- Experimental + Simulation\EM Models\Newest Model- 150209\150209-nff';
    %outputname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Actin in Infected Erythrocytes\ExtractedSurface-Different Outputs\output.txt';
    
    objects = importTriangles(fname);
    %objects = getActinRings({});
    %objects = getCrossingLines({});
    %objects = importEMData(fname);
    %objects = importPly(fname);
    %objects = importFilamentousStructures(fname);
    %objects = objects{1};
    %objects = swapColumns(objects,2,3);    
    %objects = rescaleObjects(objects,2.3); %important for mitochondria
    %objects = importSelfMeassuredMicrotubuli();
   
    [ap,ep,stormPoints,idxClust] = doSimulation(objects,bspsnm,pabs,loa,aoa,noc,doa,bspnm,rof,abpf,sxy,sz,bd,fpab,diaClust,enableClust);
    if (size(stormPoints,1)>0)
        writeOutput(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    end
    visualizeResults(objects,ap,ep,stormPoints,idxClust,cEM,cAB)
profile viewer
end

function writeOutput(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab)
    writeStormPointsForVisp(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    writeOutputFileMalk(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
    writeStormPointsForAmira(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm,fpab);
end

function [ap,ep,stormPoints,idxClust] = doSimulation(objects,bspsnm,pabs,loa,aoa,noc,doa,bspnm,rof,abpf,sxy,sz,bd,fpab,diaClust,enableClust)
    if isSurfaceData(objects)
        [ap,ep,idxClust] = findAntibodiesTri_3(objects, bspsnm, pabs, loa, aoa, noc, doa, diaClust, enableClust);
    else
        [ap,ep,~]=findLines(objects, bspnm, pabs, aoa, loa, rof);
    end
    %stormPoints = zeros(5);
    [stormPoints, ~ , ~] = findStormPoints(ep, abpf, sxy, sz, bd, fpab, true);
end

function visualizeResults(objects,ap,ep,stormPoints,idxClust,cEM,cAB)
    figHandle = findobj('Type','figure');
    close(figHandle);
    
    showEM(objects,cEM)
    showAntibodies(ap,ep,cAB)
    showEMAntibodies(objects,ap,ep,cEM,cAB)
    showClusterTriangles(objects,idxClust,cEM)
%     showStormPoints(stormPoints)
%     showEMAntibodiesStormPoints(objects,ap,ep,stormPoints,cEM,cAB)
    setViewForFigures(150,30)
end

function showClusterTriangles(objects,idxClust,col)
    triang = getMatrix(objects);
    triang = triang(:,1:3,:);
    sizeClust = size(idxClust,1);
    objectsClust = zeros(sizeClust,3,3);
    for i = 1 : sizeClust
        objectsClust(i,:,:) = triang(idxClust(i),:,:);
    end

    fig = figure(findNextFigure);
    clf
    hold on
    if isSurfaceData(objects)
        figure(fig)
        for i = 1:size(objectsClust,1)
            plot3(objectsClust(i,:,1),objectsClust(i,:,2),objectsClust(i,:,3),'color',col)
        end
    else
        printEMLines(objectsClust,fig,col)
    end
end

function showEMAntibodiesStormPoints(objects,ap,ep,stormPoints,col1,col2)
    fig = figure(findNextFigure);
    clf
    hold on
    if isSurfaceData(objects)
        printTriangles(objects,fig,col1)
    else
        printEMLines(objects,fig,col1)
    end
    printAntibodies(ap,ep,fig,col2)
    printSTORMPoints(stormPoints,fig)
end

function showStormPoints(stormPoints)
    fig = figure(findNextFigure);
    clf
    hold on
    printSTORMPoints(stormPoints,fig)
end

function showEMAntibodies(objects,ap,ep,col1,col2)
    fig = figure(findNextFigure);
    clf
    hold on
    if isSurfaceData(objects)
        printTriangles(objects,fig,col1)
    else
        printEMLines(objects,fig,col1)
    end
    printAntibodies(ap,ep,fig,col2)
end

function showAntibodies(ap,ep,col)
    fig = figure(findNextFigure);
    clf
    hold on
    printAntibodies(ap,ep,fig,col)
end

function showEM(objects,col)
    fig = figure(findNextFigure);
    clf
    hold on
    if isSurfaceData(objects)
        printTriangles(objects,fig,col)
    else
        printEMLines(objects,fig,col)
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
        try
            idx = max(figureHandles)+1;
        catch
            figNames=[figureHandles.Number];
            idx = max(figNames)+1;
        end
    end
end