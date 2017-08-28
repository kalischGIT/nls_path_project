function tot = trap_meth(dx,qvals)

    tot = dx/2*(qvals(1) + qvals(end) + 2*sum(qvals(2:end-1)));