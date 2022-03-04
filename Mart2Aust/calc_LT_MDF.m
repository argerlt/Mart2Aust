function [LT_MDF,psi] = calc_LT_MDF(CS_HT,CS_LT,halfwidth,OR,ori_count)
% create a misorientation distribution function (MDF) for the LT_Phase.
% An MDF is roughly equivilant to a probability distribution function. If
% you know the orientation relationship for the martensitic transformation
% as well as the kernel width (aka, the "spread") expected around each
% measurement (both calculated during AutoOR), you can make a theoretical
% ODF of how a single super large Austenite(HT) grain might look after
% completely transforming into Martensite(LT) variants.
% This goes one step further: if you know what the ODF looks like, you can
% ALSO make an MDF. This would let you input a misorientation (ie, the
% rotation between two orientations), and it would tell you the likelyhood
% of randomly observing that misorientation between two pixels from the
% same prior austenite grain. Higher means likely same parent, lower means
% likely different. This is what we will use to weight the out-of-plane
% weights in the graph cut later. 

%% Startup stuff. check for inputs, and switch to defaults where needed
if exist('halfwidth','var')
    assert(length(halfwidth) ==1 && isfloat(halfwidth),...
        'Invalid halfwidth. Must be a positive float');
else
    halfwidth = 0.03;
    disp('!! ==== NO HALFWIDTH PROVIDED ==== !!')
    disp('This is allowed but NOT reccomended. using default:' +halfwidth);
end

if exist('OR','var')
    assert(isvector(OR) && length(OR) ==3 && isfloat(OR),...
        'Invalid OR. Must be a vector of 3 floats')
else
    OR = [2.16,8.06,8.30]; %GT
    disp('!! ==== NO ORIENTATION RELATIONSHIP PROVIDED ==== !!')
    disp('This is allowed but NOT reccomended. using default:' +OR);
end
if ~exist('ori_count','var')
    ori_count = 1000;
    disp('calcMDF: using 1000 orientations to make MDF')
end

% Also, lets suppress some warnings cause I'm tired of seeing them
orig_warning_state = warning;
warning('off','all');
%% Create a unimodal ODF for the PAG orientations and sample it
PAG_ori = orientation(idquaternion,CS_HT);
PAG_ODF = unimodalODF(PAG_ori,'halfwidth',halfwidth);
PAG_oris = calcOrientations(PAG_ODF,ori_count);

%% For each ori, get the 24 Martensite variants and compile
R2T = calc_R2T(OR,CS_HT,CS_LT);
Mart_oris = symmetrise(PAG_oris)*R2T;
% List too big to run calcDensity on, so downsample back to ori_count
Mart_oris = Mart_oris(randperm(length(Mart_oris),ori_count));

%% Use this list and a scrambled copy to make a misorientations list and MDF
Scrambled_Mart_oris = Mart_oris(randperm(length(Mart_oris)));
% This is a quick and dirty way to get misorientation angles between
% orientations in two lists of identical lengths. read as "apply orientation B 
% to the inverse of orientation A" 
mori=inv(Mart_oris)*Scrambled_Mart_oris;
psi = deLaValleePoussinKernel('halfwidth',halfwidth);
LT_MDF=calcMDF(mori,'kernel',psi,'silent');

%Turn warnings back on
warning(orig_warning_state);

end

