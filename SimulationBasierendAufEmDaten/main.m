function main()
% profile on
    figHandle = findobj('Type','figure');
    close(figHandle);
% 	rng(5);
    rng('shuffle');
    loa = 10; %length of both antibodies combined
    aoa = 90/180*pi; %angle of antibody
    bspnm = 0.27; %binding sites per nm
    pabs = 0.40; %part of available binding sites
    abpf = 14;% average blinking per fluorophor
    rof = 11;%radius of filament

    cEM = [0,0,0]; %color of EM data
    cSTORM = [1,0,0]; %color of STORM result
    cAB = [0,1,0]; %color of Antibody
    sxy = 8; %sigma of fitting error in xy direction
    sz = 35; %sigma of fitting error in z direction

    bspsnm = .0159/2; %binding sites per square nanometer
    
    enableClust = 1;
    noc = 0.0001; %number of clusters per square nm
    areaClust = 4900; %average area of a cluster
    doa = 0.005; %denstiy of antibodies in clusters per square nm
    
    %fname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Organelle Library\Microtubules\Microtubules.wimp';
    %fname = '/media/Dev_d/Persönlicher Ordner/Docs/Skripte/Master/Studium/S_01/Kuner/EM Tomography Model/Mitochondria-Tomogram-beta-islet-cells';
    fname = '/media/Dev_d/Persönlicher Ordner/Docs/Skripte/Master/Studium/S_01/Kuner/EM Tomography Model/Mitochondria-Tomogram-beta-islet-cells.nff';
    %outputname = 'Y:\Users_shared\Superresolution Simulation Software Project- Frank and Varun\Organelle Library\Mitochondria\STORM Simulation\Mitochondria-Tomogram-beta-islet-cells.nff';
    outputname = '/media/Dev_d/Persönlicher Ordner/Docs/Skripte/Master/Studium/S_01/Kuner/EM Tomography Model/Output/Mito-Tomo_clust1_noc0.0001_areaC4900_doa0.005.nff';
    
    objects = importTriangles(fname);
    %objects = getActinRings({});
    %objects = getCrossingLines({});
    %objects = importEMData('141208-STORMmodel');
    %objects = importFilamentousStructures(fname);
    %objects = swapColumns(objects,2,3);
    %objects = rescaleObjects(objects,10);
    
    if isSurfaceData(objects)
        [ap,ep,idxClust] = findAntibodiesTri_3(objects, bspsnm, pabs, loa, aoa, noc, doa, areaClust, enableClust);
    else
        [ap,ep,idx]=findLines(objects, bspnm, pabs, aoa, loa, rof);
    end
    [stormPoints, idxF ,idxSt] = findStormPoints(ep, abpf, sxy, sz, false);

    writeStormPointsForVisp(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm);
    writeOutputFileMalk(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm);
    writeStormPointsForAmira(stormPoints,outputname,loa,aoa,bspnm,pabs,abpf,rof,sxy,sz,bspsnm);
    
    showEM(objects,cEM)
    showAntibodies(ap,ep,cAB)
    showEMAntibodies(objects,ap,ep,cEM,cAB)
    showClusterTriangles(objects,idxClust,cEM)
    showStormPoints(stormPoints)
    showEMAntibodiesStormPoints(objects,ap,ep,stormPoints,cEM,cAB)
    setViewForFigures(150,30)
% profile viewer
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
        if size(objects{i},1)>4
            issd = 0;
            break;
        end
    end
end

function idx = findNextFigure()
    figureHandles = findobj('Type','figure');
    if size(figureHandles) == 0
        idx = 1;
    else
        figNames=[figureHandles.Number];
        idx = max(figNames)+1;
    end
end
