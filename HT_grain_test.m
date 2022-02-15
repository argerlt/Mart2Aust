CS_HT = ebsd.CSList{1}; % equivilant to CS_R
CS_LT = ebsd.CSList{2}; % equivilant to CS_T
SS = ebsd.opt.SS;

searchable_ebsd = ebsd(ebsd.phaseId == 2);
for iteration = 1:10
    % find indices of a possible grain
    [Aus_Grain_Ids,grain_remainder_Ids,AusOr] = HT_grain_guess(searchable_ebsd);
    % extract the expected grain and leave the remainder for future searches
    grain_ebsd = searchable_ebsd(Aus_Grain_Ids);
    searchable_ebsd(Aus_Grain_Ids) = [];
    % remove everything from this grain that isnt the LT phase (phase 2)
    grain_ebsd = grain_ebsd(grain_ebsd.phase == 2);
    % We don't need the whole grain, just a representative sample of the
    % Martensite from a single Austenite grain. Plus, a 2x larger sample
    % causes an 8x slowdown in the graph cut, so time to downsample
    martensite=grain_ebsd.orientations;
    martensite.CS = CS_LT;
    keep=randperm(length(martensite));
    if length(martensite) > options.OR_sampling_size
        martensite=martensite(keep(1:options.OR_sampling_size));
    else
    end
    % plot the martensite pole figure if desired
    if (options.OR_plot_mart == 1)
        figure; plotPDF(martensite,Miller({0,0,1},martensite.CS),'antipodal','points','all','marker','.');
    end
%    figure()
%    aaa = searchable_ebsd(searchable_ebsd.phaseId == 2);
%    plot(aaa("Martensite"),aaa("Martensite").orientations)
    figure()
    aaa = grain_ebsd(grain_ebsd.phaseId == 2);
    plot(aaa("Martensite"),aaa("Martensite").orientations)
    
%    figure()
    ori = aaa.orientations;
    ori.CS = crystalSymmetry('m-3m');
    r = vector3d.Z;
    h = inv(ori)*r;
%    plot(h.symmetrise,'fundamentalRegion')

    figure()
    plot(h,'upper')
    
    apple = 1;
end
