function [phix,phiz] = phi_eval(theta,xi,z,Llx,eta1,K,s,cg,Om,om,k0,sig,ep)
    Kvec = pi/Llx*[0:K -K+1:-1]';
    Dx = 1i*Kvec;
    KT = 2*K;
    
    % Assume eta is in frequency space.
    eta1p = ifft(eta1);
    eta1sq = eta1p.^2;
    eta1xp = ifft(Dx.*eta1);
    meta1p = (abs(eta1p)).^2;
        
    theta2 = theta.^2;
    
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
    n2 = a1*eta1sq;
        
    zhox = -2*Om*ep*k0*meta1p;
    fhox = -s*Om*eta1p + 1i*ep*cg*eta1xp;
    shox = ep*(-2*s*Om*n2+(Om-s*om)*k0*eta1sq);
    
    fhoz = 1i*Om*eta1p + cg*ep*eta1xp;
    shoz = ep*1i*(2*Om*n2 + k0*(om - s*Om)*eta1sq);
    
    % Compute phix
    
    At10 = -ep*abs(k0)*fhox.*conj(eta1p);
    At11 = fhox.*theta;
    At12 = -ep*abs(k0)*fhox.*eta1p.*theta2;
    
    Bt10 = -ep*abs(k0)*fhoz.*conj(eta1p);
    Bt11 = fhoz.*theta;
    Bt12 = -ep*abs(k0)*fhoz.*eta1p.*theta2;
    
    %At21 = ep*(-2*abs(k0)*conj(eta1p).*shox.*theta);
    At22 = shox*theta2;  
    
    %Bt21 = ep*(-2*abs(k0)*conj(eta1p).*shoz.*theta);
    Bt22 = shoz*theta2;  
        
    Rk0 = exp(ep*z*abs(Kvec)).*fft(.5*zhox + At10);
    Rk1 = exp(z*abs(k0 + ep*Kvec)).*fft(At11);
    Rk2 = exp(z*abs(2*k0 + ep*Kvec)).*fft(At12+At22);
    phixk = exp(1i*Kvec.*xi).*(Rk0 + Rk1 + Rk2);
    
    Rk0t = exp(ep*z*abs(Kvec)).*fft(Bt10);
    Rk1t = exp(z*abs(k0 + ep*Kvec)).*fft(Bt11);
    Rk2t = exp(z*abs(2*k0 + ep*Kvec)).*fft(Bt12+Bt22);
    phizk = exp(1i*Kvec.*xi).*(Rk0t + Rk1t + Rk2t);
    
    phix = 2/KT*real( sum(phixk(1:K)) + sum(phixk(K+2:KT)) + .5*phixk(K+1) );
    phiz = 2/KT*real( sum(phizk(1:K)) + sum(phizk(K+2:KT)) + .5*phizk(K+1) );
     
end

