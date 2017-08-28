function [Tvals,xtrack,ztrack,w,sdrift] = nls_solver_pwave_1d(K,Llx,A,cg,anl,k0,Om,om,sig,ep,tf)
     
    X = (-Llx:Llx/K:Llx-Llx/K)'; % Spatial mesh
    s = sign(k0);
    n0 = om*(2*Om*s-om)/(1 + cg*om);
    disp(om*n0)
    [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
    Lind = K;
    Rind = K+2;
      
    xint1 = X(Lind);
    xint2 = X(Rind);
        
    xpos = [xint1; xint2];
    
    Omt = Om + anl*A^2*ep^2;    
    velx = @(theta) 1 + ep*k0*(om - s*Om)*(2*A*cos(theta))./(Omt + ep^2*k0*(om*n0-2*Om*k0)*A^2) ...
          + ep^2*k0*((om-2*s*Om)*a1- abs(k0)*om)*(2*A^2*cos(2*theta))./(Omt + ep^2*k0*(om*n0-2*Om*k0)*A^2);
    tvals = linspace(0,2*pi,1e3);  
    plot(tvals,velx(tvals))
    pause
    
    % Solve NLS equation in time
    opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
    [Tvals,xtrack] = ode45(@(t,lhs) phi_eval_pwave_ho_1d(t,lhs,s,A,anl,Om,om,n0,k0,sig,ep),[0 tf],xpos(2),opts);
    %[Tvals,theta_track] = ode45(@(t,lhs) phi_eval_pwave_ho_1d_theta(t,lhs,s,A,anl,Om,om,n0,a1,k0,ep),[0 tf],k0*xpos(2),opts);
    %xtrack = (theta_track - (Om + anl*A^2*ep^2)*Tvals)/k0;
    
    plot(Tvals,xtrack)
    pause
    
    sdrift = zeros(length(Tvals),1);
    ztrack = zeros(length(Tvals),1);
    
    for jj=1:length(Tvals)
        tloc = Tvals(jj);
        xloc = xtrack(jj);
        w = A*exp(1i*anl*ep^2*A^2*tloc);    
        n2 = a1*w^2;
        theta = exp(1i*(k0*xloc+Om*tloc));
        
        ztrack(jj) = ep^2*n0*A^2 + 2*ep*real( w*theta + ep*n2*theta^2);
        sdrift(jj) = Stokes_Drift_pwave(w,sign(k0),Om,om,k0);                             
    end
    
    w = A*exp(1i*anl*ep^2*A^2*tf);
    n2 = a1*w.^2;
    w = ep^2*n0*A^2 + 2*ep*real( w.*exp(1i*(k0*X + Om*tf)) + ep*n2.*exp(2*1i*(k0*X+ Om*tf)) );             
    
end