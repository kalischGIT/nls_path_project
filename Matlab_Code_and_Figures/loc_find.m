function indc = loc_find(xpos,X)
    indsl = X <= xpos;
    indsr = X > xpos;
    [~,indc] = max(indsl.*[indsr(2:end);indsr(1)]);