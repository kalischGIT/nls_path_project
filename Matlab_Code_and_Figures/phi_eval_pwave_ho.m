function rhs = phi_eval_pwave_ho(t,lhs,s,A,anl,Om,om,k0,sig,ep)
    
    x = lhs(1);
    z = lhs(2);
    
    % Assume eta is in physical space.
    theta = exp(1i*(k0*x + Om*t));
    eta1 = A*exp(1i*A^2*anl*ep^2*t);
    
    eta1sq = eta1.^2;
    eta1cu = eta1sq.*eta1;
    meta1 = A^2;
       
    q1 = -s*Om*eta1;
        
    theta2 = theta.^2;
    theta3 = theta2.*theta;
    
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
    n2 = a1*eta1sq;
    q2 = a2*eta1sq;
        
    n3 = a3*eta1cu;
    q3 = b3*eta1cu;
        
    zhox = -2*Om*ep*k0*meta1;
    fhox = q1 - ep^2*( k0^2*(om-s*Om)*eta1sq + 4*Om*k0*n2).*conj(eta1);
    shox = ep*(q2+Om*k0*eta1sq);
    thox = ep^2*(k0^2*(om-s*Om)*eta1cu + 4*Om*k0*eta1.*n2 + q3);
    
    fhoz = 1i*Om*eta1 + ep^2*1i*( -3*Om*k0^2*eta1sq + k0*om*n2 - 2*abs(k0)*Om*n2 - k0*q2 ).*conj(eta1);
    shoz = ep*1i*(2*Om*n2 + k0*(om - s*Om)*eta1sq);
    thoz = 1i*ep^2*(-3*Om*k0^2*eta1cu + 4*Om*k0*eta1.*n2 + 3*Om*n3);
        
    % Compute phix
    
    At10 = -ep*abs(k0)*fhox.*conj(eta1) + ep^2*(k0^2*fhox.*meta1);
    At11 = fhox.*theta + ep^2*(-.5*(k0^2*fhox.*eta1sq + 2*k0^2*fhox.*meta1 + 2*abs(k0)*conj(fhox).*n2).*theta);
    At12 = -ep*abs(k0)*fhox.*eta1.*theta2 + ep^2*(k0^2*fhox.*eta1sq).*theta2;
    At13 = ep^2*(-.5*fhox.*(k0^2*eta1sq + 2*abs(k0)*n2).*theta3);
    
    Bt10 = -ep*abs(k0)*fhoz.*conj(eta1) + ep^2*(k0^2*fhoz.*meta1);
    Bt11 = fhoz.*theta + ep^2*(-.5*(k0^2*fhoz.*eta1sq + 2*k0^2*fhoz.*meta1 + 2*abs(k0)*conj(fhoz).*n2).*theta);
    Bt12 = -ep*abs(k0)*fhoz.*eta1.*theta2 + ep^2*k0^2*fhoz.*eta1sq.*theta2;
    Bt13 = ep^2*(-.5*fhoz.*(k0^2*eta1sq + 2*abs(k0)*n2).*theta3);
    
    At21 = ep*(-2*abs(k0)*conj(eta1).*shox.*theta);
    At22 = shox*theta2;  
    At23 = ep*(-2*abs(k0)*eta1.*shox.*theta3);
    
    Bt21 = -2*abs(k0)*ep*conj(eta1).*shoz.*theta;
    Bt22 = shoz*theta2;  
    Bt23 = -2*abs(k0)*ep*eta1.*shoz.*theta3;
        
    Rk0 = .5*zhox + At10;
    Rk1 = exp(abs(k0)*z).*(At11+At21);
    Rk2 = exp(2*abs(k0)*z).*(At12+At22);
    Rk3 = exp(3*abs(k0)*z).*(thox+At13+At23);
    phix = 2*real(Rk0 + Rk1 + Rk2 + Rk3);
    
    Rk0t = Bt10;
    Rk1t = exp(abs(k0)*z).*(Bt11+Bt21);
    Rk2t = exp(2*abs(k0)*z).*(Bt12+Bt22);
    Rk3t = exp(3*abs(k0)*z).*(thoz+Bt13+Bt23);
    phiz = 2*real(Rk0t + Rk1t + Rk2t + Rk3t);
       
    xdot = om*z + ep*phix;
    zdot = ep*phiz;
    rhs = [xdot;zdot];
end

