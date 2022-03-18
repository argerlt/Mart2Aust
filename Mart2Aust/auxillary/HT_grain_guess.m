function [Aus_Grain_Ids,grain_remainder_Ids,AusOr] = HT_grain_guess(ebsd,OR_guess)
% combination of MTEX-driven and graph cut `reconstruction' algorithm used
% to segment a portion of a grain expected to be from a single Austenite
% grain, which can then be used for the automatic determination of the OR

ksiKS = [5.26,10.30,10.53];
ksiNW = [0,9.74,9.74];
CS_LT = ebsd.CSList{2};
CS_HT = ebsd.CSList{1};

% Allow the use of an initial Orientation Relationship guess. Otherwise,
% default to GT (ask eric, not sure what this is
if exist('OR_guess','var')
    assert(isvector(OR_guess) && length(OR_guess) ==3 && isfloat(OR_guess),...
        "OR_guess must be a length 3 array of floats")
    ksi = OR_guess; % Begin with provided guess
else
    ksi = [2.16,8.06,8.30];     % Begin with GT measurement
end


% Compute variants and corresponding groupoid from euler angles, and
% convert the groupoids to Euler angles
[V,~] = YardleyVariants(ksi);
[G,~] = GroupoidVariantsCubic(V);
Ax = vector3d(G(:,2),G(:,3),G(:,4));
G_eul = orientation('axis',Ax,'angle',G(:,1),CS_LT);

% Classify groupoid as a misorientation in MTEX terms
rot = orientation('Euler',[0,0,0],CS_LT);
G_misos(1) = inv(rot) * rot;
G_misos(2:length(G_eul)+1) = inv(rot) * G_eul;
G_misos = transpose(G_misos);

% Set tolerance (higher for KS and NW since expirimental usuallly varys
% from expected in those cases)
if (all(ksi==ksiKS) || all(ksi==ksiNW))
    tol = 1.25*degree;
else
    tol = 0.85*degree;
end

% make a copy of ebsd we can goof with
tmpEbsd = ebsd(ebsd.phaseId == 2);
tmpEbsd.phaseMap = [1,2,0];
tmpEbsd.CSList = {tmpEbsd.CSList{1} tmpEbsd.CSList{2} tmpEbsd.CSList{end}};

% use the MTEX calcgrains to segment out Martensite grains
[grains,tmpEbsd.grainId] = calcGrains(tmpEbsd);

% Isolate the grain boundaries and adjust the phase to singular
gB = grains.boundary;
gB.phase= 2;

% Find the austenite grain boundaries based on the above tolerance
gB_aus = [];
for i = 1:length(G_misos)
    if i == 1
        gB_aus = gB(angle(gB.misorientation,G_misos(i)) < tol);
    else
        gB_aus = vertcat(gB_aus,gB(angle(gB.misorientation,G_misos(i)) < tol));
    end
end
% Merge grains together based on the chosen austenite grain boundaries
[mergedGrains,parentID]=merge(grains,gB_aus,'calcMeanorientation');
mergedGrains.CS = CS_LT;

% Find the maximum `austenite' grain and extract the Ebsd coordinates
% related to the original dataset
[~,maxGrnid] = max(mergedGrains.grainSize);
ParGrnId = find(parentID==maxGrnid);
EbId = [];
for i = 1:length(ParGrnId)
    EbId = vertcat(EbId,find(tmpEbsd.grainId==ParGrnId(i)));
end


% Truncated ebsd dataset from `grain reconstruction'
tmpEbsd = tmpEbsd(EbId);
tmpEbsd.phase = 2;

% Compute misorientation distribution function based on KS-OR and
% corresponding groupoid
hw = 2*degree;
psi=deLaValleePoussinKernel('halfwidth',hw);
%MartModf = calcODF(G_misos,'kernel',psi);
MartModf = calcDensity(G_misos,'kernel',psi);

% Establish transformation matrices
[T2R,~]=calc_T2R(ksi,CS_HT,CS_LT);
[R2T,~]=calc_R2T(ksi,CS_HT,CS_LT);

% From our EBSD region, transform the martensite to austenite and find the
% modal austenite orientation to use as our guess
%aus_odf=calcODF(symmetrise(tmpEbsd.orientations)*T2R,'kernel',psi);
aus_odf=calcDensity(symmetrise(tmpEbsd.orientations)*T2R,'kernel',psi);
[~,AusOr] = max(aus_odf);

% Call function to compute adjacency array and corresponding
% misorientations for the neighbgoring points
[adjpts,mori,~] = adjpt_moris(tmpEbsd,1);

% In-plane and out-of-plane weights
modf_vals=eval(MartModf,mori);
modf_vals(modf_vals<0)=0;

% If the input OR is KS or NW, adjust the parameterization for these. If
% it's GT or something similar, alter as well and include for twins.
if (all(ksi==ksiKS) || all(ksi==ksiNW))
    IP = [10,15];
    OP = [2,1.5e-1];
    %    TwnWts = zeros(length(tmpEbsd),1);
else
    IP = [7,8];
    OP = [2,1.5e-1];
    %    Twns = vector3d([1,1,1;-1,-1,1;-1,1,1;1,-1,1]');
    %    Twn_eul = orientation('axis',Twns,'angle',60*degree,ebsd{2}.CS);
    %    AusOrRot = rotation('euler',AusOr.phi1,AusOr.Phi,AusOr.phi2);
    %    TwnAus = AusOrRot * Twn_eul;
    %    TwnOrs = symmetrise(TwnAus)*R2T;
    %    TwnODF = calcODF(TwnOrs,'kernel',psi);
    %    TwnWts = (eval(TwnODF,tmpEbsd.orientations)+OP(1)).*OP(2);
end

% Compute ODF for parent and establish the second set of out-of-plane
% weights
From_OrGuess=symmetrise(AusOr)*R2T;
From_GuessODF=calcODF(From_OrGuess,'kernel',psi);
OP2wts=(eval(From_GuessODF,tmpEbsd.orientations)+OP(1)).*OP(2);

% Establish constant as first set of out-of-plane weights
% OP1wts = abs(1./OP2wts);
% OP1wts(OP1wts==Inf)=1e2;
% OP1wts = OP1wts + TwnWts;

IPwts = (modf_vals+IP(1)).*IP(2);

% Add source node
DiGraph=digraph;

% Establish constant as first set of out-of-plane weights
OP1wts = abs(4./OP2wts);
OP1wts(OP1wts==Inf)=1e2;

% Call function to set up the initial graph
[DiGraph,endnode,sinknode] = graph_setup(DiGraph,adjpts,IPwts,OP1wts,2);

% Add edges corresponding to the guess reconstruction orientation
DiGraph=rmedge(DiGraph,(1:length(tmpEbsd))+endnode,sinknode);
DiGraph=addedge(DiGraph,(1:length(tmpEbsd))+endnode,sinknode,OP2wts);

% Perform graph cut
[~,~,~,ct]=maxflow(DiGraph,1,sinknode);

% Truncate ct and finalize the cut
ct(end) = [];
Cut = ct - endnode;

% Convert the IDs in the grain to the Ids of the original scan and return
Aus_Grain_Ids = EbId(Cut);
% EbAus = tmpEbsd(Cut);
grain_remainder_Ids = EbId;
grain_remainder_Ids(Cut) = [];
end