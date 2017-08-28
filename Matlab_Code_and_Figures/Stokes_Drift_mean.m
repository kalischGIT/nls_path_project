function ud = Stokes_Drift_mean(xi,Llx,eta1,K,s,Om,om,k0)

    Kvec = pi/Llx*[0:K -K+1:-1]';
    KT = 2*K;
        
    % Assume eta is in frequency space.
    
    q1 = -s*Om*eta1;
    qt = 1i*Om*eta1;
    
    % Compute phix
    R1 = q1;
    intx = exp(1i*Kvec*xi).*R1;  
    phix = (sum(intx(1:K)) + sum(intx(K+2:end)) + .5*intx(K+1))/KT;  
    
    % Compute phiz    
    R1t = qt;
    intz = exp(1i*Kvec*xi).*R1t;
    phiz = (sum(intz(1:K)) + sum(intz(K+2:end)) + .5*intz(K+1))/KT;    
    
    % Compute Stokes Drift Terms
        
    ud = -2*k0/Om*(abs(phix)^2 + abs(phiz)^2 + om/Om*imag(conj(phix).*phiz));
       
end

