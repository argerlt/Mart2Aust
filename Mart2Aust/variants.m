function int_map = variants(orig_ebsd,HT_ebsd,options,LT_MDF)
% split the grains into individual variant sections.
% Process

% Pre-flight stuff:
% If users provide a misorientation distribution function for the low-temp
% phase (LT_MDF), overwrite the saved one with that
if exist('LT_MDF','var')
    orig_ebsd.opt.LT_MDF = LT_MDF;
end
Orientation_relationship = orig_ebsd.opt.OR;
psi = orig_ebsd.opt.psi;
LT_CS = HT_ebsd.CSList{2};

% Make a copy of the original that ONLY includes successfully reconstructed
% areas (ie, areas where the HT_ebsd phase is phase 3)
LT_ebsd = orig_ebsd(HT_ebsd.phaseId == 3);
assert(all(LT_ebsd.phaseId == 2)) % double check the scans line up.
% Now mask out the HT_ebsd the same way, so both contain EXACTLY identical
% data (assuming HT_ebsd was generated from a reconstruction of orig_ebsd)
HT_ebsd = HT_ebsd(HT_ebsd.phaseId == 3);
% rotate each reconstructed HT grain by the parent orientation, thus
% creating a map of the varient orientations relative to their parent.
v_ori = inv(rotation(HT_ebsd.orientations)).* LT_ebsd.orientations;


% Calculate all the data for the OP weights
% NOTE: for this next bit, I make 24 Orientation Distribution functions,
% and evaluate v_ori for every one. this could be done other ways (rotate
% var map, single eval, etc) but this seems to work fine.
[V,~] = YardleyVariants(Orientation_relationship);
Var_weights = zeros(size(v_ori,1),size(V,1));
for i = 1:size(V,1)
    V_eul = orientation('matrix',transpose(V{i}),LT_CS);
    %        Var_odf = calcDensity(symmetrise(V_eul),'kernel',psi);
    Var_odf = calcDensity(V_eul,'kernel',psi);
    Var_weights(:,i) = eval(Var_odf,v_ori); %#ok<EV2IN> 
end
Var_weights(Var_weights<=0) = 0;
% reweight them using the y=mx+1 with the appropriate weights from options
Var_weights = (Var_weights.*options.Seg_OP_m)+options.Seg_OP_b;

% Create MODF for orientations expected within a given block. Use this to
% Calculate the IP weights (Ask steve/eric about why groupoids are used)
[Gtmp,~] = GroupoidVariantsCubic(V);
G0 = orientation('euler',[0,0,0],LT_CS);
Ax = vector3d(Gtmp(:,2),Gtmp(:,3),Gtmp(:,4));
G_eul = [G0;orientation('axis',Ax,'angle',Gtmp(:,1),LT_CS)];
clear Gtmp G0 Ax
% Compute MODF for variant boundaries
variant_modf = calcDensity(symmetrise(G_eul(1:2),LT_CS),'kernel',psi);
neighborhood_size = int32(options.degree_of_connections_for_neighborhood);
neigh_list = find_neighboring_pairs(LT_ebsd,neighborhood_size);
IP_wts = get_In_plane_weights(neigh_list,LT_ebsd,options,variant_modf);

% create list of unique Aus orientations, along with a maskable grain map
[unique_parent_oris,~,grain_map] =  HT_ebsd.orientations.unique;
unique_grain_count = size(unique_parent_oris,1);
% make an array for storing the variant ids
int_map = (1:size(LT_ebsd,1))*0;

% now loop through all the grains
for i = 1:unique_grain_count
    grain_mask = grain_map == i;
    LT_act = LT_ebsd(grain_mask);
    VW_act = Var_weights(grain_mask,:);
    [NL_act,IP_act] = prune_IP_graph_connections(LT_act.id,neigh_list,IP_wts);
    int_map_act = find_grain_variants(LT_act,VW_act,NL_act,IP_act);
    int_map(grain_mask) = int_map_act +(i*24);
    fprintf('grain %d of %d\n',i,unique_grain_count)
end
untransformed_count = sum(int_map == 0);
if untransformed_count >0
    disp('=========================================')
    disp('Warning: Segmentation did not complete!!!')
    fprintf('%0.0f voxels remain untransformed\n',untransformed_count)
    disp('=========================================')
else
    disp('=========================================')
    disp('Segmentation completed successfully!')
    disp('=========================================')
end

end

%% ====================================================================
function    int_map = find_grain_variants(LT_act,VW_act,NL_act,IP_act)
% use min cut/max flow to make per-pixel variant orientation assignments.
N = size(LT_act,1);
L = NL_act(:,1);
R = NL_act(:,2);
assigned = 10^(ceil(log10(max(max(VW_act))))+4);
for ii = 1:2:22
    in_rows = ii:ii+1;
    out_rows  = setdiff(1:size(VW_act,2),in_rows);
    OP_in_block = max(VW_act(:,in_rows),[],2);
    OP_not_in_block = max(VW_act(:,out_rows),[],2);
    % Denoise (maybe delete this later? not convinced it is good to do)
    % The "Denoise" refers to making every pixel the average of its 4
    % nearest neighbors to smooth out erratic values in the OP weights
    % Commenting out because slow
    %OP_in_block_DN = OP_in_block*0;
    %OP_not_in_block_DN = OP_not_in_block * 0;
    %for i = 1:size(LT_act,1)
    %    OP_in_block_DN(i) = mean(OP_in_block([R(L == i);L(R == i)]));
    %    OP_not_in_block_DN(i) = mean(OP_not_in_block([R(L == i);L(R == i)]));
    %end

    % make a digraph with n+2 nodes (1 per voxel,plus source and sink)
    % NOTE: source has ID n+1, sink has ID n+2
    FP_digraph = digraph;
    FP_digraph = addnode(FP_digraph,N+2);

    % add source-to-voxel weights (equal to likelyhood that pixel DOES belong
    % to a given block)
    FP_digraph = addedge(FP_digraph, N+1, 1:N, OP_in_block);
    % add voxel-to-sink weights (equal to likelyhood that a pixel belongs to
    % one of the blocks not yet cut out)
    FP_digraph = addedge(FP_digraph, 1:N, N+2, OP_not_in_block);

    % Add in-plane (voxel to voxel) connections)
    FP_digraph = addedge(FP_digraph,L,R,IP_act);
    FP_digraph = addedge(FP_digraph,R,L,IP_act);

    % Perform graph cut
    [~,~,cs,~]=maxflow(FP_digraph,N+1,N+2);
    cs(cs>length(LT_act)) = [];

    if numel(cs)>0
    % for the parts IN the block, decide which variant each pixel is by
    % figuring out which has the highest likelyhood
%    [~,q] = min(abs(VW_act(cs,in_rows)-OP_in_block(1,cs)'),[],2);
    [~,q] = min(abs(VW_act(cs,in_rows)-OP_in_block(cs,1)),[],2);
    % set all weights for the extracted parts to zero
    VW_act(cs,:) = 0;
    % assign the variants an extremely high OP weight
    VW_act(cs(q ==1),ii) = assigned;
    VW_act(cs(q ==2),ii+1) = assigned;
    % cut all IP weights attached to those points so they don't get recut
    IP_act(ismember(L,cs)) = 0;
    IP_act(ismember(R,cs)) = 0;
    end
end

% for the last block, just choose the max of the remaining two variants.
% This makes sure EVERY pixel gets assigned and avoids some edge cases that
% can cause breaks
[~,q] = max(VW_act,[],2);
VW_act(q>22,:) = 0;
VW_act(q==23,23) = assigned;
VW_act(q==24,24) = assigned;

% make the int map to pass back to the segmentation routine
[~,int_map] = max(VW_act,[],2);
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


function IP_wts = get_In_plane_weights(neigh_list,ebsd,options,MODF)
% use the misorientation between neighbor pairs to determine a "weight" to
% give to that connection in terms of network flow capacity.

% get misorientation angle between pairs
[~,id_Dl] = ismember(neigh_list(:,1),ebsd.id);
[~,id_Dr] = ismember(neigh_list(:,2),ebsd.id);
o_Dl = ebsd(id_Dl).orientations;
o_Dr = ebsd(id_Dr).orientations;

%o_Dl = ebsd(neigh_list(:,1)).orientations;
%o_Dr = ebsd(neigh_list(:,2)).orientations;
Mori = inv(o_Dl).*(o_Dr);

% Find likelyhoods for those misorientation angles to occur
Mori.SS = MODF.SS;
MDF_vals=eval(MODF,Mori); %#ok<EV2IN> 
MDF_vals(MDF_vals<0)=0;

% Alter their weights using a y=mx+b style linear equation. For anyone
% reading this and confused, we are altering the "strength" of the
% pixel-to-pixel (ie, in-plane) connections. increasing the value of b
% increases the cost of making any in-plane cuts (ie, will favor a shorter
% overall length of grain boundaries). increasing m will increase the
% cost of making a cut through a likely pair compared to an unlikely pair.
m = options.Seg_IP_m;
b = options.Seg_IP_b;
IP_wts = MDF_vals*m +b;
end
