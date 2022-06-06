function [packet, block, variant] = SGID_to_PBV(subgrain_id)
%UNTITLED3 takes in a list of subgrain IDs and returns the relavent 
% Variant, Block, and Packet ID.

% make variant ID
variant = rem(subgrain_id,24)+1;
% make block ID map
block = floor((variant+1)/2);
% make packed ID map
packet = floor((variant+5)/2);

end