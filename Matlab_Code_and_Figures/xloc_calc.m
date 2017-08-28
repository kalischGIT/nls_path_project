function xloc = xloc_calc(xpos,cgsc,t,Tval)

    xloctemp = Tval/pi*(xpos + cgsc*t) - Tval;
    
    if xloctemp > Tval
        xloc = -Tval + (xloctemp-Tval);
    elseif xloctemp < -Tval
        xloc = Tval - (-xloctemp-Tval);
    else
        xloc = xloctemp;
    end