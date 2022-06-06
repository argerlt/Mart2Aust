% ======================================================================= %
%       Plotting.m
%       Shows some basic plotting methods
% ======================================================================= %
%% Check if a reconstruction was already ran, and if not, run "Basic.m"
if size(dir('Original.ang'),1) ~= 1 ||...
        size(dir('Reconstruction.ang'),1) ~= 1 ||...
        size(dir('AF96.txt'),1) ~= 1
    Basics
end

%% Load reconstruction data
% Load original ebsd
O = EBSD.load('Original.ang');
% Load reconstructed ebsd
R = EBSD.load('Reconstruction.ang');
% load reconstruction and plotting options
options = load_options('Steel');
% Load variant map
likelihood = table2array(readtable('AF96.txt','Range','G:G'));
sg_map = table2array(readtable('AF96.txt','Range','H:H'));
[p,b,v] = SGID_to_PBV(sg_map);
% Load the OR and HW from the "Original.ang" scan (there HAS to be a better
% way to do this...)
temp = split(string(fileread('Original.ang')),"# ");
OR = split(temp(startsWith(temp,'OR')));
OR = transpose(str2double(OR(2:4)));
HW = split(temp(startsWith(temp,'HW')));
HW = transpose(str2double(HW(2)));
O.opt.OR = OR;
O.opt.HW = HW;

% MTEX tries to intelligently rewrite the PhaseID order when it load ebsd
% files. This backfires for us, since we used a custom PhaseID ordering.
% The correct fix is to rewrite prep_for_recon, Auto_OR, and subgrain to
% use an ordering MTEX automatically understands. However, the way it
% orders changed between 5.1.1 and 5.7.0, and will probably change again.
% the quick and easy fix is to just use a different phase ID ordering in
% the loaded ebsd scans:
%
% 1 : Unindexed
% 2 : Untransformed Parent (High temperature, or HT)
% 3 : Transformed Child (Low Temperature, or LR)
% 4 : Reconstructed Parent (starts empty) (Reconstructed, or R)
% 5 : Reconstructed subgrain (starts empty) (subgrain and/or variant, V)
Mart = O(O.phaseId == 3);
[Mart_grains,Mart.grainId] = calcGrains(Mart,'angle',5*degree);
PAG = R(R.phaseId == 4);
[PAG_grains,PAG.grainId] = calcGrains(PAG,'angle',2*degree);

%% Make maps of the grain and subgrain boundaries. 
% NOTE:FOR THIS TO LOOK CORRECT, WE ARE GOING FILL AND SMOOTH SCANS.
% THEREFORE,THESE EBSD OBJECTS SHOULD ONLY BE USED FOR VISUALIZATION, NOT
% CALCULATIONS.
Martf = O(Mart(Mart_grains(Mart_grains.grainSize >= 5)).id);
[Martf_grains,Martf.grainId] = calcGrains(Martf,'angle',5*degree);

PAGf = R(PAG(PAG_grains(PAG_grains.grainSize >= 10)).id);
[PAGf_grains,PAGf.grainId] = calcGrains(PAGf,'angle',2*degree);

% Its also handy to have a subset of the gB that are in the twinning 
% condition
gB_PAG = PAGf_grains.boundary('Reconstructed Austenite','Reconstructed Austenite');
% S9 done from Physical Metallurgy
CS = PAG.CS;
u1 = Miller( 1, 2, 1,CS);
v1 = Miller(-1,-2,-1,CS);
u2 = Miller(-1, 0, 1,CS);
v2 = Miller( 1, 0,-1,CS);
S3_twin_miso = orientation('map',u1,v1,u2,v2);
% S9 was done from 1984 Hans Grimmer Paper
S9_angle = 38.94*degree;
S9_axis = vector3d(1,1,0);
S9_twin_miso = orientation('axis',S9_axis,'angle',S9_angle);

S3_isTwinning = angle(gB_PAG.misorientation,S3_twin_miso) < 0.5*degree;
S9_isTwinning = angle(gB_PAG.misorientation,S9_twin_miso) < 0.5*degree;
S3_TwinBoundaries = gB_PAG(S3_isTwinning);
S9_TwinBoundaries = gB_PAG(S9_isTwinning);
All_Twin_Boundaries = [S9_TwinBoundaries S3_TwinBoundaries];

%% Packet, Block, and Variant coloring
key = ipfHSVKey(Mart);
mp = key.orientation2color(Mart.orientations);
steel_cmap = get_subgrain_colormaps('Steel');

%% Reconstruction Plots
% Plot the original scan
figure()
plot(Mart, Mart.orientations)
% same plot but no micron bar
figure()
plot(Mart,Mart.orientations,'micronbar','off')

% Plot the Reconstructed microstructure
figure()
plot(PAG, PAG.orientations,'micronbar','off')

% Plot the original microstructure but with the PAG outlined
figure()
plot(Mart,Mart.orientations,'micronbar','off')
hold on
plot(PAGf_grains.boundary,'linecolor','k','linewidth',1.5)

% Plot the original microstructure but with the PAG outlined and twinning
% boundaries highlighted
figure()
plot(Mart,Mart.orientations,'micronbar','off')
hold on
plot(gB_PAG,'linecolor','k','linewidth',6,'displayName','Prior')
plot(S3_TwinBoundaries,'linecolor','red','linewidth',2,'displayName','S3 Twin')
plot(S9_TwinBoundaries,'linecolor','blue','linewidth',2,'displayName','S9 Twin')
hold off

% Plot the Reconstructed microstructure but with the PAG outlined and twinning
% boundaries highlighted
figure()
plot(PAGf,PAGf.orientations,'micronbar','off')
hold on
plot(gB_PAG,'linecolor','k','linewidth',4,'displayName','Prior')
plot(S3_TwinBoundaries,'linecolor','red','linewidth',1.5,'displayName','S3 Twin')
plot(S9_TwinBoundaries,'linecolor','blue','linewidth',1.5,'displayName','S9 Twin')
hold off

%% Plot Austenite Likelihood Plot

% Plot likelihood plot of chosen austenite orientations transforming into 
% the observed martensitic microstructure
figure()
plot(O,likelihood,'MicronBar','off')

% likelihood plot but with the PAG outlined. (Useful for highlighting bad
% reconstructions)
figure()
plot(O,likelihood,'MicronBar','off')
hold on
plot(gB_PAG,'linecolor','k','linewidth',3,'displayName','Prior')
plot(S3_TwinBoundaries,'linecolor','red','linewidth',1.5,'displayName','S3 Twin')
plot(S9_TwinBoundaries,'linecolor','blue','linewidth',1.5,'displayName','S9 Twin')
hold off
%% Plot Packets

%% Plot Blocks

%% Plot Variants
figure()
for i = 1:24
    mask = rem(v(Mart.id),24) ==(rem(i+1,24));
    %mask = logical(mask.*transpose(Mart.grainId ==mode(Mart.grainId)));
    plot(Mart(mask),mp(mask),'FaceColor',steel_cmap.Variant(i,:))
    hold on
end

%% Plot single grain with packets and Blocks outlined

% Plot grain structure or merged twins with the corresponding numbering of
% each grain (if set to 1; otherwise, it will just plot the grain boundary
% structure without labeling
%  plotGrains(myEBSD.AusGrains,1);
