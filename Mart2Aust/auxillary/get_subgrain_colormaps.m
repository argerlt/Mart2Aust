function cmap = get_subgrain_colormaps(system)
% currently only works for steel. since Aus/Mart has hierarchical coloring,
% returns colormapping for Variants, Blocks, and Packets
if system == 'Steel'
    steel_cmap.Packet = [
        1,    0,    0;
        0,    1,    0;
        0,    0,    1;
        1,    1,    0];
    steel_cmap.Block = [
        252,  122,  120 ;   250,  4,    0;
        102,  0,    0   ;   120,  252,  122;
        4,    250,  0   ;   0,    125,  2;
        130,  125,  253 ;   2,    2,    250;
        1,    1,    100 ;   220,  220,  170;
        250,  250,  2   ;   227,  176,  10
        ]./255;

    steel_cmap.Variant = [
        227,    97,    95;   255,   147,   145;
        225,     4,     0;   255,     4,     0;
        77,     0,     0;   127,     0,     0;
        95,   227,    97;   145,   255,   147;
        4,   225,     0;     4,   255,     0;
        0,   100,     2;     0,   150,     2;
        105,   100,   228;   155,   150,   255;
        2,     2,   225;     2,     2,   255;
        1,     1,    75;     1,     1,   125;
        195,   195,   145;   245,   245,   195;
        225,   225,     2;   255,   255,     2;
        202,   151,    10;   251,   201,    10
        ]./255;

    cmap = steel_cmap;
    return
else
    error('currently, only steel is implimented.')
end