filename = "Data//4D-XIII-A_41out.ang";
filename = "../EBSD/AF96_321x/AF_001.ang";
ebsd = EBSD.load(filename,'convertEuler2SpatialReferenceFrame');

%ebsd = EBSD.load(join(split(filename,"/"),filesep),'convertEuler2SpatialReferenceFrame');
%ebsd = EBSD.load(filename,'convertEuler2SpatialReferenceFrame');
plotx2east

ebsd.CSList{2}.mineral = 'Austenite'
ebsd.CSList{3}.mineral = 'Martensite'

%% Load it up and group into initial grains
% grain reconstruction
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 3*degree);
% remove small grains
ebsd(grains(grains.grainSize < 3)) = [];
% reidentify grains with small grains removed:
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',3*degree);
grains = smooth(grains,5);
% plot the data and the grain boundaries
figure(1); 
plot(ebsd('Martensite'),ebsd('Martensite').orientations,'figSize','large')
figure(2); 
plot(ebsd('Martensite'),ebsd('Martensite').orientations,'figSize','large')
hold on
plot(grains.boundary,'linewidth',2)
hold off

%% set up the recon job and use it to optimize the OR

job = parentGrainReconstructor(ebsd,grains);

% initial guess for the parent to child orientation relationship
job.p2c = orientation.KurdjumovSachs(job.csParent, job.csChild)
job.calcParent2Child
figure(3);
histogram(job.fit./degree)
xlabel('disorientation angle')

%% Compute grain-to-grain misorientation, and use that to seed adjacency matrix in recon "job"
% compute the misfit for all child to child grain neighbours
[fit,c2cPairs] = job.calcGBFit;
% select grain boundary segments by grain ids
[gB,pairId] = job.grains.boundary.selectByGrainId(c2cPairs);

% plot the child phase, and on top of it the boundaries colorized by the misfit
figure(4)
plot(ebsd('Martensite'),ebsd('Martensite').orientations,'figSize','large','faceAlpha',0.5)
hold on;
% scale fit between 0 and 1 - required for edgeAlpha
plot(gB, 'edgeAlpha', (fit(pairId) ./ degree - 2.5)./2 ,'linewidth',2);
hold off


%% actual reconstruction line here
job.calcGraph('threshold',2.5*degree,'tolerance',2.5*degree);

figure(5)
plot(ebsd('Martensite'),ebsd('Martensite').orientations,'figSize','large','faceAlpha',0.5)
hold on;
job.plotGraph('linewidth',2)
hold off

% cluster the grains in the graph into logical parent grains
job.clusterGraph('inflationPower',2.6)

% compute parent orientations
job.calcParentFromGraph

% plot them
figure(6)
plot(job.parentGrains,job.parentGrains.meanOrientation)

% plot misfit as well
figure(7)
plot(job.grains,job.grains.fit./degree)
%plot(job.grains, job.grains.clusterSize < 15)
setColorRange([0,5])
mtexColorbar

% NOTE: theere is another step here inthe example i am not fully
% comprehending, where the user can revert child grains with large misfits
% and then add them to different parent grains, which could be useful for
% us, but I doubt it
%% POST (This is where this matters again)

% merge grains with similar orientation
job.mergeSimilar('threshold',7.5*degree);

% plot the result
figure(8)
plot(job.parentGrains,job.parentGrains.meanOrientation)

job.calcVariants

% associate to each packet id a color and plot
color = ind2color(job.transformedGrains.packetId);
figure(9)
plot(job.transformedGrains,color,'faceAlpha',0.5)
hold on
parentGrains = smooth(job.parentGrains,2);
plot(parentGrains.boundary,'linewidth',3)

% outline a specific parent grain
grainSelected = parentGrains(parentGrains.findByLocation([100,80]));

hold on
plot(grainSelected.boundary,'linewidth',3,'lineColor','w')
hold off


% identify childs of the selected parent grain
childGrains = job.grainsPrior(job.mergeId == grainSelected.id);

% plot these childs
figure(10)
plot(childGrains,childGrains.meanOrientation)

% and top the parent grain boundary
hold on
plot(grainSelected.boundary,'linewidth',2)
hold off



% the measured child orientations that belong to parent grain 279
childOri = job.ebsdPrior(childGrains).orientations;

% the orientation of parent grain 279
parentOri = grainSelected.meanOrientation;

% lets compute the variant and packeIds
[variantId, packetId] = calcVariantId(parentOri,childOri,job.p2c);

% colorize child orientations by packetId
color = ind2color(packetId);
figure(11)
plotPDF(childOri,color, Miller(0,0,1,childOri.CS),'MarkerSize',2,'all')

% the positions of the parent (001) directions
hold on
plot(parentOri.symmetrise * Miller(0,0,1,parentOri.CS),'markerSize',10,...
  'marker','s','markerFaceColor','w','MarkerEdgeColor','k','linewidth',2)

% the theoretical child variants
childVariants = variants(job.p2c, parentOri);
plotPDF(childVariants, 'markerFaceColor','none','linewidth',1.5,'markerEdgeColor','k')
hold off



parentEBSD = job.calcParentEBSD;

% plot the result
plot(parentEBSD('Austenite'),parentEBSD('Austenite').orientations,'figSize','large')

