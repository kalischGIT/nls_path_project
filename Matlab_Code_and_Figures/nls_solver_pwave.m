function [Tvals,xtrack,ztrack,zdiff,w,sdrift] = nls_solver_pwave(K,Llx,A,anl,k0,Om,om,sig,ep,tf)
     
    X = (-Llx:Llx/K:Llx-Llx/K)'; % Spatial mesh
    s = sign(k0);
    
    uint = A;
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
    n2 = a1*uint^2;
    n3 = a3*uint^3;
    
    sint = 2*ep*real( uint*exp(1i*k0*X) + ep*n2*exp(2*1i*k0*X) + ep^2*n3*exp(3*1i*k0*X));
                    
    Lind = K;
    Rind = K+2;
      
    xint1 = X(Lind);
    xint2 = X(Rind);
        
    xpos = [xint1; xint2];
    zpos = [sint(Lind);sint(Rind)];
                      
    % Solve NLS equation in time
    opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
    [Tvals,Y2] = ode45(@(t,lhs) phi_eval_pwave_ho(t,lhs,s,A,anl,Om,om,k0,sig,ep),[0 tf],[xpos(2);zpos(2)],opts);
    xtrack = Y2(:,1);
    ztrack = Y2(:,2);
    sdrift = zeros(length(Tvals),1);
    zdiff = zeros(length(Tvals),1);
    
    for jj=1:length(Tvals)
        tloc = Tvals(jj);
        w = A*exp(1i*anl*ep^2*A^2*tloc);    
        
        sdrift(jj) = log10(abs(Stokes_Drift_pwave(ztrack(jj),w,sign(k0),Om,om,k0,ep))); 
        
        n2 = a1*w^2;
        n3 = a3*w^3;
        wp = 2*ep*real( w*exp(1i*(k0*X+ Om*tloc)) + ep*n2*exp(2*1i*(k0*X+ Om*tloc)) + ep^2*n3*exp(3*1i*(k0*X+ Om*tloc)) );
            
        yy = spline(X,wp,xtrack(jj));
        zdiff(jj) = log10(abs(ztrack(jj)-yy));              
    end
    
    w = A*exp(1i*anl*ep^2*A^2*tf);
    n2 = a1*w.^2;
    n3 = a3*w.^3;
    w = 2*ep*real( w.*exp(1i*(k0*X + Om*tf)) + ep*n2.*exp(2*1i*(k0*X+ Om*tf)) + ep^2*n3.*exp(3*1i*(k0*X+ Om*tf)) );        
        
end