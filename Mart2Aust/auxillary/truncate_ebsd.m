function [trunc_ebsd] = truncate_ebsd(ebsd,xmin,xmax,ymin,ymax)

    rr = [xmin ymin (xmax-xmin) (ymax-ymin)];
    cond = inpolygon(ebsd,rr);
    trunc_ebsd = ebsd(cond);

end