function neighboring_pairs = find_neighboring_pairs(ebsd,connection_distance)
% given a list of X, Y locations, and a unit cell, this function
% calculates an adjacency matrix.

% Austin hotfix: this is SUPER slow on big (1000x1000) ebsd scans, AND it
% scales geometrically. so, we split this up into slices instead.
unique_x = unique(sort(unique(ebsd.x)));

all_neighboring_pairs = [];
for i  = 1:100:size(unique_x,1)
    l_mask = ebsd.x >= unique_x(i);
    h_mask = ebsd.x <= (unique_x(min(i+104,size(unique_x,1))));
    mask = l_mask.*h_mask.*[1:size(ebsd,1)]';
    mask = mask(mask>0);
    sliced_ebsd = ebsd(mask);
    % Note: I(Austin) don't fully understand this function. I_FD is called 'the
    % incidence matrix between faces to cells', but I can confirm that matrix
    % multiplying it DOES give you a valid adjacency matrix (ie if you want to
    % know if voxel A touches voxel B, go to the coordinates (A,B). If the
    % value there is non-zero, then they touch
    [~,~,I_FD] = spatialDecomposition([sliced_ebsd.prop.x,sliced_ebsd.prop.y],sliced_ebsd.unitCell);
    A = I_FD.' * I_FD;

    if ~exist('connection_distance','var')
        connection_distance = 1;
    end
    assert(connection_distance == 1 ||connection_distance == 2 ||connection_distance == 3,...
        "The only allowable neigh_count values are 1 (only 1st nearest neighbors,"+...
        "2 (neighbors of neighbors), or 3 (neighbors of neighbors of neighbors)")
    if connection_distance > 1
        A = A +A*A;
    end
    if connection_distance > 2
        A = A +A*A;
    end

    % extract just the point pairs where A_D != 0 (ie, coordinates the
    % correspond to a pair of neighbors
    [Dl, Dr] = find(A);
    % only document each set of pairs once, and remove self connections
    use = Dl > Dr;
    Dl = Dl(use);
    Dr = Dr(use);

    % At this point, Dl and Dr are correct neighbor pairs, but the Ids they
    % refer to are the indicies of the voxel in the ebsd obhect, NOT the voxels
    % id. We need to do a quick retranslation to fix this

    ids = sliced_ebsd.id;
    Dl_trans = ids(Dl);
    Dr_trans = ids(Dr);

    % repackage
    neighboring_pairs = [Dl_trans,Dr_trans];
    all_neighboring_pairs = [all_neighboring_pairs ;neighboring_pairs];
end
neighboring_pairs = unique(all_neighboring_pairs,'rows');
end

