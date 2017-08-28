function [phix,phiz] = phi_eval_pwave(theta,A,z,n1,s,Om,om,k0,sig,ep)

    q1 = -s*Om*n1;
     
    theta2 = theta.^2;
    
    q2 = 2*Om*k0*(2*Om^2-2*sign(k0)*om*Om+om^2)/(4*Om^2-2*sign(k0)*(k0*(1+4*sig*k0^2)+om*Om));
    n2 = abs(k0)*(2*Om^2-2*sign(k0)*om*Om+om^2)/(4*Om^2-2*sign(k0)*(k0*(1+4*sig*k0^2)+om*Om));
       
    qt = 1i*Om*n1;
    
    % Compute phix
    R1 = theta*q1.*exp(z*abs(k0));
    R20 = -Om*k0*A^2/pi;
    
    R21 = theta2*exp(2*z*abs(k0))*(q2+k0*Om)*n1^2;
    
    R22 = -s*Om*abs(k0)*(theta2*exp(2*z*abs(k0))*n1^2  + A^2);

    phix = 2*real( R1 + ep*(R21 - R22) ) + ep*R20;  
    
    % Compute phiz    
    R1t = theta*qt*exp(z*abs(k0));
    
    R20 = 1i*(k0*(s*Om-om) + Om*abs(k0))*A^2;
    
    R22 = 1i*theta2*(2*Om*n2-k0*(s*Om-om))*n1^2*exp(2*z*abs(k0));
    
    R2ps = -1i*Om*theta2*n1^2*exp(2*z*abs(k0));
          
    tvec = R1t + ep*(R20+R22+R2ps);
    
    phiz = 2*real(tvec);    
        
end

