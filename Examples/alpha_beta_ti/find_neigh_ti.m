function [final_neighboring_pairs] = find_neigh_ti(ebsd, phase_id)
% FIND_NEIGH_TI given a list of X, Y locations, and a unit cell, this 
% function calculates an adjacency matrix.
%   Detailed explanation goes here
unique_x = unique(sort(unique(ebsd.x)));

all_neighboring_pairs = [];
% spatial decomp is slow, so lets split the math up into overlapping chunks
% of 100 lines
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

    % extract just the point pairs where A_D != 0 (ie, coordinates the
    % correspond to a pair of neighbors
    [Dl, Dr] = find(A);
    % only document each set of pairs once, and remove self connections
    use = Dl > Dr;
    Dl = Dl(use);
    Dr = Dr(use);

    % repackage
    neighboring_pairs = [Dl,Dr];
    all_neighboring_pairs = [all_neighboring_pairs ;neighboring_pairs];
end
% At this point, Dl and Dr are correct neighbor pairs, but the Ids they
% refer to are the indicies of the voxel in the ebsd obhect, NOT the voxels
% id. We need to do a quick retranslation to fix this

single_phase_ebsd = ebsd(ebsd.phase == phase_id);
ids = single_phase_ebsd.id;
Dl_trans = ids(all_neighboring_pairs(:,1));
Dr_trans = ids(all_neighboring_pairs(:,2));

final_neighboring_pairs = [Dl_trans,Dr_trans];

end