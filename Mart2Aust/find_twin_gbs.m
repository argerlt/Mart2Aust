function [Aus_gb_merged,S3,S9] = find_twin_gbs(ebsd, phaseID)
% Returns all the grain boundaries, as well as the S3 and S9 boundaries

% find all boundaries
[grains,~,~] = calcGrains(ebsd, 'angle', 1*degree);
grains = grains.smooth(5);
Recon_phasename = ebsd.CSList{3}.mineral;
Aus_gb = grains.boundary(Recon_phasename,Recon_phasename);

% define S3 and S9. 
%S3 is defined using u1v1 method from Steve
CS_A = ebsd.CSList{3};
u1 = Miller( 1, 2, 1,CS_A);
v1 = Miller(-1,-2,-1,CS_A);
u2 = Miller(-1, 0, 1,CS_A);
v2 = Miller( 1, 0,-1,CS_A);
S3_twin_miso = orientation('map',u1,v1,u2,v2);
% S9 was done from 1984 Hans Grimmer Paper
S9_angle = 38.94*degree;
S9_axis = vector3d(1,1,0);
S9_twin_miso = orientation('axis',S9_axis,'angle',S9_angle,CS_A);
S9_twin_miso.SS = S9_twin_miso.CS;

% restrict to twinnings with threshold 0.5 degree
S3_isTwinning = angle(Aus_gb.misorientation,S3_twin_miso) < 1.5*degree;
S9_isTwinning = angle(Aus_gb.misorientation,S9_twin_miso) < 1.5*degree;

S3 = Aus_gb(S3_isTwinning);
S9 = Aus_gb(S9_isTwinning);
All_Twin_Boundaries = [S9 S3];

[Aus_gb_merged,~] = merge(grains,All_Twin_Boundaries);

end