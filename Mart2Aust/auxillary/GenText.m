function [] = GenText(name, original, recon, likelihoods, sub_grains)
%UNTITLED2 Creates M2A text output from an original reconstructed EBSD map
% plus a variant map.
%"Name" should be the text string for the filename being written
% "original" should be the input ebsd
% "recon" should be the reconstructed ebsd
% "likelihood" should be the likelihood map from Reconstruction
% "sub_grains" should be the sub_grains map. will calculate if not given
if nargin < 4
    fprintf("no variant map given, calculating instead")
    sub_grains = variants(sub_grains);
end


% take the subgrain map and add in all the non-reconstructed points as -1
sg_map = zeros(size(recon,1),1) -1;
sg_map(recon.phase == 3) = sub_grains;
% get the variant, block, and packet maps
[p, b, v] = SGID_to_PBV(sg_map);

% store all the data in a table with headers
T = table(original.rotations.phi1,...
    original.rotations.Phi,...
    original.rotations.phi2,...
    recon.rotations.phi1,...
    recon.rotations.Phi,...
    recon.rotations.phi2,...
    likelihoods,sg_map, v, b, p,...
    'VariableNames',...
    {'phi1','Phi','phi2',...
    'Recon_phi1','Recon_Phi','Recon_phi2',...
    'likelihoods','subgrain_ID','Variant','Block','Packet'});
% save the data as a text file
writetable(T,name);
end
