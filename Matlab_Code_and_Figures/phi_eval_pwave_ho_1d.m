function rhs = phi_eval_pwave_ho_1d(t,lhs,s,A,anl,Om,om,n0,k0,sig,ep)
    
    x = lhs;
    
    % Assume eta is in physical space.
    theta = exp(1i*(k0*x + Om*t));
    eta1p = A*exp(1i*A^2*anl*ep^2*t);
    
    eta1sq = eta1p.^2;
             
    theta2 = theta.^2;
    
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
    n2 = a1*eta1sq;
        
    z = ep^2*n0*A^2 + 2*ep*real( eta1p*theta + ep*n2*theta2);
    ztrun = 2*ep*real( eta1p*theta );      
    fhox = -s*Om*eta1p;
    shox = ep*(-2*s*Om*a1+(Om-s*om)*k0)*eta1sq;
   
     % Compute phix
    
    At11 = fhox.*theta;
    At12 = -ep*abs(k0)*fhox.*eta1p.*theta2;
    
    At22 = shox*theta2;  
    
    Rk1 = (1 + abs(k0)*ztrun).*At11;
    Rk2 = (At12+At22);
    phix = 2*real(Rk1 + Rk2);
    
    rhs = om*z + ep*phix;
    
end

