O = Tasks(9).ebsd;
R = Tasks(9).Recon_ebsd;
options = Tasks(9).options;

MO = O(O.phaseId == 2);
[gmap, MO.grainId] = calcGrains(R(R.phaseId == 3));

v_map = variants(O, R, options);
key = ipfHSVKey(MO);
mp = key.orientation2color(MO.orientations);
steel_cmap = get_subgrain_colormaps('Steel');

figure()
for i = 1:24
    mask = rem(v_map,24)==(rem(i,24));
    mask = logical(mask.*transpose(MO.grainId ==mode(MO.grainId)));
    plot(MO(mask),mp(mask),'FaceColor',steel_cmap.Variant(i,:))
    hold on
end