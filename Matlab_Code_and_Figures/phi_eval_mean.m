function phimx = phi_eval_mean(xi,z,Llx,eta1,K,s,Om,om,k0,ep)
    Kvec = pi/Llx*[0:K -K+1:-1]';
    KT = 2*K;
    
    % Assume eta is in frequency space.
    eta1p = ifft(eta1);
    meta1p = (abs(eta1p)).^2;
    
    cp = -Om/k0;
            
    q1 = -s*Om*eta1;
    q1p = ifft(q1);
        
    zhox = -2*Om*k0*meta1p;    
    fhox = q1p;
    fhoz = 1i*Om*eta1p;
    
    At10 = -abs(k0)*fhox.*conj(eta1p);
    At11 = fhox;
    
    Bt11 = fhoz;
    
    %Rk0 = exp(1i*Kvec.*xi).*exp(ep*z*abs(Kvec)).*fft(.5*zhox + At10);
    %Rk1 = exp(1i*Kvec.*xi).*exp(z*abs(k0 + ep*Kvec)).*fft(At11);
    
    %Rk1t = exp(1i*Kvec.*xi).*exp(z*abs(k0 + ep*Kvec)).*fft(Bt11);
    
    Rk0 = exp(1i*Kvec.*xi).*fft(.5*zhox + At10);
    Rk1 = exp(1i*Kvec.*xi).*fft(At11);
    
    Rk1t = exp(1i*Kvec.*xi).*fft(Bt11);
        
    R0 = 2/KT*real( sum(Rk0(1:K)) + sum(Rk0(K+2:KT)) + .5*Rk0(K+1) );
    R1 = 1/KT*( sum(Rk1(1:K)) + sum(Rk1(K+2:KT)) + .5*Rk1(K+1) );
    
    R1t = 1/KT*( sum(Rk1t(1:K)) + sum(Rk1t(K+2:KT)) + .5*Rk1t(K+1) );
    
    phimx = R0 + 2*( real( R1*conj(R1) + R1t*conj(R1t) ) + om/(k0*(cp-om*z))*imag(R1*conj(R1t)) )/(cp-om*z);
   
end

