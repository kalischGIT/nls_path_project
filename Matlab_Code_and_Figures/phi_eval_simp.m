function [phix,phiz] = phi_eval_simp(theta,Llx,indf,eta1,pfac,K,s,cg,Om,om,k0,sig,ep)
    Kvec = pi/Llx*[0:K -K+1:-1]';
    Dx = 1i*Kvec;
    
    % Assume eta is in frequency space.
    eta1p = ifft(eta1);
    eta1sq = eta1p.^2;
    eta1cu = eta1sq.*eta1p;
    eta1xp = ifft(Dx.*eta1);
    eta1peta1px = eta1p.*eta1xp;
        
    q1 = -s*Om*eta1;
    q1p = ifft(q1);
        
    theta2 = theta.^2;
    theta3 = theta2.*theta;
    
    a1 = abs(k0)*(2*Om^2-2*s*om*Om+om^2)/(4*Om^2-2*s*(k0*(1+4*sig*k0^2)+om*Om));
    a2 = 2*Om*k0*(2*Om^2-2*s*om*Om+om^2)/(4*Om^2-2*s*(k0*(1+4*sig*k0^2)+om*Om));
    
    n2 = a1*eta1sq;
    q2 = a2*eta1sq;
        
    a3t = -1/2*(3*k0^5*sig - 2*Om^2*k0^2*s - 2*Om*a1*abs(k0)*om + 4*Om^2*a1*k0 - 2*Om*a2*abs(k0) + 2*Om*k0^2*om + 2*a1*k0*om^2 + 2*a2*k0*om);
    b3t = -(2*Om*a1*abs(k0) - k0^2*s*om + a2*k0);
    
    n3 = (s*a3t/3 - Om*b3t)/(s*(9*k0^3*sig+k0+om*Om)/3-Om^2)*eta1cu;
    q3 = (-Om*a3t+(9*k0^3*sig+k0+om*Om)*b3t)/(s*(9*k0^3*sig+k0+om*Om)/3-Om^2)*eta1cu;
    
    q0 = -(om-2*s*Om)*real(ifft(abs(Dx).*fft(abs(eta1p).^2)));
    n0 = -cg*q0./(1+cg*(1+om));        
        
    etaf = ep^2*n0 + 2*real(eta1p.*theta + ep*n2.*theta2 + ep^2*n3.*theta3);        
    Qf = ep^2*q0 + 2*real(q1p.*theta + ep*q2.*theta2 + ep^2*q3.*theta3);   
    
    etafx = 2*real(1i*k0*eta1p.*theta + 2*1i*k0*ep*n2.*theta2 + 3*1i*k0*ep^2*n3.*theta3);
    etaft = 2*real(1i*Om*eta1p.*theta + 2*1i*Om*ep*n2.*theta2 + 3*1i*Om*ep^2*n3.*theta3);
    
    etafX = 2*real(ep*eta1xp.*theta);
    etafT = 2*real(ep*cg*eta1xp.*theta + 2*a1*cg*ep^2*eta1peta1px.*theta2);
    etaftau = 2*real(ep^2*pfac*eta1xp.*theta);
    
    phix = Qf(indf) - ep*(etafx(indf) + ep*etafX(indf))*(etaft(indf) + ep*etafT(indf)) ...
           - ep^2*(om*etaf(indf) + Qf(indf))*(etafx(indf))^2;
    
    phiz = etaft(indf) + ep*etafT(indf) + ep^2*etaftau(indf) - ep^2*(etafx(indf))^2*etaft(indf)...
           + ep*(Qf(indf) + om*etaf(indf))*(etafx(indf) + ep*etafX(indf));
