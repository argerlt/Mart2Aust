function [Aus_gb_merged,ebsd,S3,S9] = find_twin_gbs(ebsd, phaseID)
    % Returns all the grain boundaries, as well as the S3 and S9 boundaries,
    % and stores the grain id in the ebsd data set
    
    % find all boundaries
    [grains,ebsd.grainId,~] = calcGrains(ebsd,'unitCell','angle', 1*degree);
    grains = grains.smooth(5);
    Recon_phasename = ebsd.CSList{3}.mineral;
    Aus_gb = grains.boundary(Recon_phasename,Recon_phasename);
    % NOTE 1: Alex add in: 'unitCell' doesn't allow for grains to
    % eat up other grains (phases, unindexed points, etc)
    
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
    
    % Merge parents and twins
    [Aus_gb_merged,ParentId] = merge(grains,All_Twin_Boundaries);
    
    % Loop through ebsd and replace grainIds with merged ones
    for ii = 1:length(grains.id)
        ebsd.grainId(ismember(ebsd.grainId,ii)) = ParentId(ii);
    end
    
    % Assign parent orientations to the merged grain structure.
    for jj = 1:length(Aus_gb_merged.id)
        MergedGrns  = grains.id(ismember(ParentId,jj));
        [~,MaxId]  = max(grains.grainSize(MergedGrns));
        ParOr       = grains(MergedGrns(MaxId)).meanOrientation;
        Aus_gb_merged(jj).meanOrientation = ParOr;
    end
    % NOTE 2: When you merged the grains, a lot of the grain orientations 
    % were assigned values of nan. This fixes that, but under the 
    % assumption that the largest surface area for each identified parent 
    % and twin is the parent, and all others are the twins. This isn't 
    % necessarily correct, but can be fixed later on and is fast

end