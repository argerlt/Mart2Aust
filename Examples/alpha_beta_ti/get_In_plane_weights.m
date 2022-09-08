function IP_wts = get_In_plane_weights(neigh_list,ebsd,MDF,m,b)
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

% Find likelyhoods for those misorientation angles to occur in the expected
% Martensite structure
LT_MDF_vals=eval(MDF,Mori);
LT_MDF_vals(LT_MDF_vals<0)=0;

% Alter their weights using a y=mx+b style linear equation. For anyone
% reading this and confused, we are altering the "strength" of the
% pixel-to-pixel (ie, in-plane) connections. increasing the value of b
% increases the cost of making any in-plane cuts (ie, will favor a shorter
% overall length of grain boundaries). increasing m will increase the
% cost of making a cut through a likely pair compared to an unlikely pair.
IP_wts = LT_MDF_vals*m +b;
%look up orientation of IDs in neighborhood list, find the angle between
% them, look up the likelyhood for that weight in the misorientation
% distribution function, and Finall

end
