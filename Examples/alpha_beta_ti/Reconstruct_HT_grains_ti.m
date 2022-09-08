function [LT_ebsd, likelihoods] = Reconstruct_HT_grains_ti(LT_ebsd,...
    neigh_list, IP_wts)
%Reconstruct_HT_grains perform the actual reconstruction
%   performs a 'while ' loop that progressively cuts out possible grains
%   until eithr there is nothing left to cut, or the algorithm starts
%   failing too much
Active_Ebsd = LT_ebsd(LT_ebsd.phase == 2);
continue_recon = true;
iterations = 0;
bad_cut_counter = 0;
likelihoods = zeros(Active_Ebsd.size);
%temp plotting stuff
%close all;
%%figure()
%%plot(Active_Ebsd,Active_Ebsd.orientations)

while continue_recon == 1
    % checks to break out of while loop
    if iterations > LT_ebsd.opt.max_recon_attempts || bad_cut_counter > 10
        continue_recon = false;
        disp('\n=========================================')
        disp('    reconstruction failed to complete    ')
        fprintf('    %0.0f voxels remain untransformed\n',sum(LT_ebsd.phase == 2))
        disp('=========================================\n')
        continue
    end
    if sum(LT_ebsd.phase == 2) == 0
        continue_recon = false;
        disp('\n=========================================')
        disp('reconstruction completed successfully!')
        disp('=========================================\n')
        continue
    end

    %% ======== Step 1 ======== %%
    % Find a likely high temp grain orientation to attempt to cut out.
    Guess_ID = randsample(1:length(Active_Ebsd),1,true);
    HT_guess_ori = Active_Ebsd(Guess_ID).orientations;

    HT_guess_ori.CS = Active_Ebsd.CSList{2};

    % Using that guess, make a rough cut to find a likely HT grain.
    Rough_Guess = First_Pass_Cut(Active_Ebsd,neigh_list,IP_wts,HT_guess_ori);
    % is the grain big enough? if not, iterate counters and try again.
    if size(Rough_Guess,1) < LT_ebsd.opt.min_cut_size
        bad_cut_counter = bad_cut_counter +1;
        iterations = iterations+1;
        continue
    end

    %% ======== Step 2 ======== %%
    % Now we want to clean up this guess with a precision cut, where we
    % let the code choose the most likely parent orientation instead of
    % guessing
    [proposed_grain, Guess_ori] = Precision_cut(Rough_Guess,neigh_list,IP_wts);
    % is the grain big enough? if not, iterate counters and try again.
    if size(proposed_grain,1) < LT_ebsd.opt.min_cut_size
        bad_cut_counter = bad_cut_counter +1;
        iterations = iterations+1;
        continue
    end

    % if the grain made it this far, reset the counters.
    bad_cut_counter = 0;
    LT_ebsd(proposed_grain.id).orientations = proposed_grain.orientations;
    LT_ebsd(proposed_grain.id).phase = 12;
    % Also assign the likelihoods relative to each twin to the likelihood
    % map
    likelihoods(proposed_grain.id) = eval(...
        LT_ebsd.opt.parent_child_MDF,proposed_grain.orientations);
    Active_Ebsd = LT_ebsd(LT_ebsd.phase == 2);

    % report on reconstruction progress
    iterations = iterations +1;
    fprintf('\n ------ Iter: %d Pcnt: %0.2f Remainder:%0.0f ------\n',...
        iterations,...
        sum(LT_ebsd.phase == 2)*100./size(LT_ebsd,1),...
        sum(LT_ebsd.phase == 2))
end
% At this point, either the reconstruction is finished, or it failed to
% complete. either way, send the final result back to Aus_Recon

end


%% ====================================================================
function Rough_Guess = First_Pass_Cut(Active_ebsd,neigh_list,IP_wts,HT_guess_ori,options)
% Does a single layer graph cut to pull out an area that likely contains at
% least one prior austenite grain, plus some surrounding materials

% do some setup stuff
%%%HT_guess_ori = orientation.byEuler(296.5*degree,231.8*degree,49.4*degree,Active_ebsd.CS) % erase this later
%%%%%HT_guess_ori = orientation.byEuler(68*degree,140.5*degree,185.5*degree,Active_ebsd.CS); % erase this later
[pruned_neigh_list,pruned_IP_wts] = prune_IP_graph_connections(Active_ebsd.id,neigh_list,IP_wts);
N = size(Active_ebsd,1); % number of voxels
L = pruned_neigh_list(:,1); % left side of neigh_list connections
R = pruned_neigh_list(:,2);% right side of neigh_list connections
oris = Active_ebsd.orientations;
mori = inv(HT_guess_ori).*oris;
LT_MDF = Active_ebsd.opt.LT_MDF;
%clear Active_ebsd
likelyhoods = eval(LT_MDF,mori);
likelyhoods(likelyhoods <=0) = 0;
%lll = ((likelyhoods)+2)*0.175;

% temp plotting stuff
%%figure()
%%plot(Active_ebsd,lll)

% make a digraph with n+2 nodes (1 per voxel,plus source and sink)
% NOTE: source has ID n+1, sink has ID n+2
FP_digraph = digraph;
FP_digraph = addnode(FP_digraph,N+2);

% add source-to-voxel weights (equal to likelyhood that a given voxel's
% orientation came from the suggested prior austenite ODF, as calculated
% from the transformed phase's MDF)
% NOTE: eval calls take time, so we should make as few as possible. However,
% this one changes every loop. Could pregen several lists from the start,
% but that would be memory intensive so probably not worth it.
%likelyhoods = eval(LT_MDF,mori);
%likelyhoods(likelyhoods <=0) = 0;
%OP_wts = (likelyhoods*options.RGC_post_pre_m) + options.RGC_post_pre_b;
OP_wts = (likelyhoods*3) + 12;
% OP_wts = (likelyhoods +2)*0.175;
FP_digraph = addedge(FP_digraph, N+1, 1:N, OP_wts);

% Add in-plane (voxel to voxel) connections)
FP_digraph = addedge(FP_digraph,L,R,pruned_IP_wts);
FP_digraph = addedge(FP_digraph,R,L,pruned_IP_wts);

% add voxel-to-sink weights (all equal to the mean of the OP weights)
FP_digraph = addedge(FP_digraph, 1:N, N+2, ones(N,1)*mean(OP_wts));
%FP_digraph = addedge(FP_digraph, 1:N, N+2, 4./OP_wts);

% Perform graph cut
[~,~,cs,~]=maxflow(FP_digraph,N+1,N+2);
cs(cs>length(Active_ebsd)) = [];
Rough_Guess = Active_ebsd(cs);
% Code for debugging to show the cut out area. NOTE: this is not a grain,
% Its just a region of the scan that likely has at least one complete grain
% in it
%%figure();
%%plot(Rough_Guess,Rough_Guess.orientations)
% Here is some extra troubleshooting code for seeing the IP weights, which
% are what this algorithm is trying to cut along (blue = weak,
% yellow=strong)
%%figure()
%%l = pruned_neigh_list(:,1);
%%r = pruned_neigh_list(:,2);
%%scatter(-Active_ebsd(l).y -Active_ebsd(r).y, Active_ebsd(l).x +Active_ebsd(r).x,1, pruned_IP_wts)
end



function [proposed_grain, Guess_ori] = Precision_cut(Rough_Guess,neigh_list,IP_wts,options)
% Starting with a rough prior cut, this cut finds the most common high temp
% orientation and cuts out JUST that grain and and twins of it.

% do some setup stuff (note this is a more heavily pruned starting list
% than the rough cut)
[pruned_neigh_list,pruned_IP_wts] = prune_IP_graph_connections(Rough_Guess.id,neigh_list,IP_wts);
N = size(Rough_Guess,1); % number of voxels
L = pruned_neigh_list(:,1); % left side of neigh_list connections
R = pruned_neigh_list(:,2);% right side of neigh_list connections
psi = Rough_Guess.opt.pc_kernel;
oris = Rough_Guess.orientations;

% find the most likely HT orientation, deduce (if applicable), generate all
% LT variants of the parent and twin (120 non-unique for steel), and build
% an ODF from them whose kernel spread matches the estimated per-variant
% spread of the LT phase (determined in the Auto-OR script). This is how we
% will weight the out of plane weights. (STEVE: speed this up if you can, 
% huge time sink. approx 10-30% of run total time on the next few lines)
alpha2beta = Rough_Guess.opt.alpha2beta;
Possible_beta_oris = oris.symmetrise*inv(alpha2beta);
Parent_odf=calcDensity(Possible_beta_oris,'kernel',psi);
% maybe modes?
Guess_ori = calcModes(Parent_odf);
beta2alpha = Rough_Guess.opt.beta2alpha;

expected_variants =symmetrise(Guess_ori)*inv(beta2alpha);
variant_odf = calcDensity(expected_variants,'kernel',psi);

likelyhoods = eval(variant_odf,oris);
likelyhoods(likelyhoods <=0) = 0;
OP_wts = (likelyhoods*12) + 3;

% make a digraph with n+2 nodes (1 per voxel,plus source and sink)
% NOTE: source has ID n+1, sink has ID n+2
FP_digraph = digraph;
FP_digraph = addnode(FP_digraph,N+2);

% add source-to-voxel weights (equal to likelyhood that a given voxel's
% orientation is part of the grain of the suggested orientation)
FP_digraph = addedge(FP_digraph, N+1, 1:N, OP_wts);
% add voxel-to-sink weights (PRECISION CUT: WE ARE NOW WEIGHTING BY
% INVERSE OF SOURCE-TO-VOXEL WEIGHTS)
FP_digraph = addedge(FP_digraph, 1:N, N+2, 4./OP_wts);

% Add in-plane (voxel to voxel) connections)
FP_digraph = addedge(FP_digraph,L,R,pruned_IP_wts);
FP_digraph = addedge(FP_digraph,R,L,pruned_IP_wts);


% Perform graph cut
[~,~,cs,~]=maxflow(FP_digraph,N+1,N+2);
cs(cs>length(Rough_Guess)) = [];
proposed_grain = Rough_Guess(cs);
% Code for debugging to show the cut out area. NOTE: this is not a grain,
% Its just a region of the scan that likely has at least one complete grain
% in it
%%figure();
%%plot(proposed_grain,proposed_grain.orientations)
% Here is some extra troubleshooting code for seeing the IP weights, which
% are what this algorithm is trying to cut along (blue = weak,
% yellow=strong)
%figure()
%l = pruned_neigh_list(:,1);
%r = pruned_neigh_list(:,2);
%scatter(-Active_ebsd(l).y -Active_ebsd(r).y, Active_ebsd(l).x +Active_ebsd(r).x,1, IP_wts)
end


function [Parent_and_Twin_IDs,likelihoods] = Seperate_Twins(proposed_grain, neigh_list, IP_wts, PT_oris,options)
% Given we have a for sure parent grain, check to see if parts would make
% more sense as a twin or as part of the parent.

% do some setup stuff
[pruned_neigh_list,pruned_IP_wts] = prune_IP_graph_connections(proposed_grain.id,neigh_list,IP_wts);
N = size(proposed_grain,1); % number of voxels
L = pruned_neigh_list(:,1); % left side of neigh_list connections
R = pruned_neigh_list(:,2);% right side of neigh_list connections
OR = proposed_grain.opt.OR;
HT_CS = proposed_grain.CS;
[R2T,~] = calc_R2T(OR,proposed_grain.CSList(3),proposed_grain.CSList(2));
psi = proposed_grain.opt.psi;
oris = proposed_grain.orientations;

% We already got the most likely HT Parent grain. NOTE: We COULD
% recalculate the most likely grain here (as was done in the original
% code), but I am fairly sure it is mostly redundant to do so.

Parent_and_Twin_IDs = zeros(size(PT_oris,1),size(oris,1));

% Find the likelyhood that each pixel is part of the Parent grain. this
% will become the background weighting for the 4 twin cuts.
Parent_variants = rotation(symmetrise(PT_oris(1)))*R2T;
Parent_odf = calcDensity(Parent_variants,'kernel',psi);
Parent_odf.CS = HT_CS;
likelihoods = eval(Parent_odf,oris);%#ok<EV2IN> 
likelihoods(likelihoods <=0) = 0;
mx = (likelihoods*options.RGC_post_pre_m);
b =  options.RGC_post_pre_b;
OP_Parent_wts = mx + b;
clear Parent_variants Parent_odf mx b

for i = 2:size(PT_oris,1)

    % find the likelyhood that each pixel is part of grain i
    Twin_variants = rotation(symmetrise(PT_oris(i)))*R2T;
    Twin_odf = calcDensity(Twin_variants,'kernel',psi);
    Twin_odf.CS = HT_CS;
    fg_likelyhoods = eval(Twin_odf,oris);%#ok<EV2IN> 
    fg_likelyhoods(fg_likelyhoods<0) = 0;
    mx = (fg_likelyhoods*options.RGC_post_pre_m);
    b =  options.RGC_post_pre_b;
    OP_twin_wts = mx + b;
    clear Twin_variants Twin_odf mx b
    
    % Set the weights of already assigned poitns to 1e-10, making already
    % assigned grains nearly impossible to cut (should maybe be 0?)
    assigned =(1:N).*(sum(Parent_and_Twin_IDs,1)>0);
    OP_twin_wts(ismember((1:N),assigned)>0) = 1e-10;

    % make a digraph with n+2 nodes (1 per voxel, plus source and sink)
    % NOTE: source has ID n+1, sink has ID n+2
    FP_digraph = digraph;
    FP_digraph = addnode(FP_digraph,N+2);
    % add source-to-voxel weights (likelyhood pixel is twin i)
    FP_digraph = addedge(FP_digraph, N+1, 1:N, OP_twin_wts);
    % add voxel-to-sink weights (likelyhood pixel is part of the parent)
    FP_digraph = addedge(FP_digraph, 1:N, N+2, OP_Parent_wts);
    % Add in-plane (voxel to voxel) connections)
    FP_digraph = addedge(FP_digraph,L,R,pruned_IP_wts);
    FP_digraph = addedge(FP_digraph,R,L,pruned_IP_wts);
    % Perform graph cut
    [~,~,cs,~]=maxflow(FP_digraph,N+1,N+2);
    cs(cs>length(proposed_grain)) = [];
    % save out the mask of the cut
%    yes_mask = [1:N].*ismember([1:N],cs);
    Parent_and_Twin_IDs(i,1:end) = ismember((1:N),cs);
    % find the already assigned orientations, set their in-plane weights to
    % zero so they won't get pulled out again.
    assigned =(1:N).*(sum(Parent_and_Twin_IDs,1)>0);
    neigh_mask = ((ismember(L,assigned) == 0) +(ismember(R,assigned) == 0))>0;
    pruned_IP_wts = pruned_IP_wts.*neigh_mask;    
    likelihoods = max(fg_likelyhoods.*transpose(Parent_and_Twin_IDs(i,1:end)),likelihoods);

    % plotting for help
    %figure()
    %plot(proposed_grain,OP_wts)
    %figure()
    %l = L(L>0);
    %r = L(R>0);
    %scatter(-proposed_grain(l).y -proposed_grain(r).y, ...
    %    proposed_grain(l).x +proposed_grain(r).x,1, pruned_IP_wts)
    %%%figure()
    %%%plot(proposed_grain,Parent_and_Twin_IDs(i,1:end))
end
% get a map of twice-twinned (should never happen) and never twinned pixels
overlap_mask = sum(Parent_and_Twin_IDs,1) <= 1;
parent_mask = sum(Parent_and_Twin_IDs,1) ~= 1;
% Twice twinned pixels get switched to never twinned
Parent_and_Twin_IDs = Parent_and_Twin_IDs.*overlap_mask;
% assign never twinned pixels to parent
Parent_and_Twin_IDs(1,:) = parent_mask;
% and add the IDs back in (maybe make this step unnecessary in the future?)
Parent_and_Twin_IDs = Parent_and_Twin_IDs.*(1:N);
end

function [remapped_neigh_list,pruned_IP_wts] = prune_IP_graph_connections(ids,neigh_list,IP_wts)
% take the global neighbor list and in plane weights, prune out the
% connections that don't connect voxels in the ids list, and renumber them
% appropriately so the digraph function doesn't make orphan nodes
ids = sort(unique(ids));
prune_mask = ismember(neigh_list(:,1),ids).*ismember(neigh_list(:,2),ids);

% prune
pruned_IP_wts = IP_wts(prune_mask == 1);
pruned_neigh_list = neigh_list(prune_mask == 1,:);


% create translater array
translator = zeros(max(ids),1);
for i = 1:length(ids)
    translator(ids(i))=i;
end

%remap
l = translator(pruned_neigh_list(:,1));
r = translator(pruned_neigh_list(:,2));

remapped_neigh_list = [l r];


end

function OP_wts = get_Out_of_plane_weights(oris,guess_ori,MDF,options)
% given a list of orientations, guess ori, and MODF, use eval to find the
% likelyhood that each orientation came from the guess orientation, then
% weight those values according to the OP weighting and scaling factors

% get misorientation angle
mori = inv(guess_ori).*oris;
% query MDF for likelyhood of those misorientation angles
likelyhoods = alt_eval(MDF,mori);
likelyhoods(likelyhoods <=0) = 0;
% apply a y = mx+b style rescaling factor. Increasing m makes close
% likelyhoods easier to include and unlikely grains harder. increasing b
% just makes ALL cuts equally harder (only really matters for precision
% cut, where some of the weights are done as 1/OP)
OP_wts = (likelyhoods*options.RGC_post_pre_m)+options.RGC_post_pre_b;


end

function orphan_IDs = prune_orphaned_pixels(Active_Ebsd,neigh_list)
% filter out pixels without any neighbors

a = neigh_list(ismember(neigh_list(1:end,1),Active_Ebsd.id),2);
b = neigh_list(ismember(neigh_list(2:end,1),Active_Ebsd.id),1);
ab = unique([a; b]);
orphan_IDs = Active_Ebsd.id(ismember(Active_Ebsd.id,ab) == 0);
end
    % BIG NOTE HERE: FOR SOME REASON, the original code uses not just
    % neighbors, but neighbors of neighbors in their code (IE,1st through
    % 3rd 2D Von Neumann neighborhoods) to calculate the adjacency matrix. 
    % This SEEMS wrong, as it isnt a typical adjacency array A, but the 
    % equivilant of A +(A*A). It would also slow down the graph cut by a 
    % factor of 4. HOWEVER, for some reason it works, so I'm keeping as an 
    % option.To switch to first (IE, only 1st Von Neumann
    % neighborhood) change options.degree_of_connections_for_neighborhood
    % from 2 to 1.


