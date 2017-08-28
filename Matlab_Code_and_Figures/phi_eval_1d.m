function phix  = phi_eval_1d(theta,xi,z,Llx,eta1,K,s,cg,Om,om,k0,sig,ep)
    Kvec = pi/Llx*[0:K -K+1:-1]';
    Dx = 1i*Kvec;
    KT = 2*K;
    
    % Assume eta is in frequency space.
    eta1p = ifft(eta1);
    eta1sq = eta1p.^2;
    eta1xp = ifft(Dx.*eta1);
        
    theta2 = theta.^2;
    
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
    n2 = a1*eta1sq;
        
    fhox = -s*Om*eta1p + 1i*ep*cg*eta1xp;
    shox = ep*(-2*s*Om*n2+(Om-s*om)*k0*eta1sq);
    
    % Compute phix
    
    At11 = fhox.*theta;
    At12 = -ep*abs(k0)*fhox.*eta1p.*theta2;
    At22 = shox*theta2;  
    
    Rk1 = exp(z*abs(k0 + ep*Kvec)).*fft(At11);
    Rk2 = exp(z*abs(2*k0 + ep*Kvec)).*fft(At12+At22);
    phixk = exp(1i*Kvec.*xi).*(Rk1 + Rk2);
    
    phix = 2/KT*real( sum(phixk(1:K)) + sum(phixk(K+2:KT)) + .5*phixk(K+1) );
     
end

