function [xtrack,ztrack,w,sdrift] = nls_stepper_1d(K,Llx,ad,anl,cg,k0,Om,om,sig,ep,tf,dt)

    % Llx is domain size, i.e. solve on domain [-Llx,Llx]
    % sig is -1 for bright/focusing, +1 for dark/defocusing
    % tf is total simulation time
    % mu is magnitude of damping on envelope of dark initial condition
        nmax = round(tf/dt); % Step size for time integrator and number of time steps
        Tval = Llx;
        X = (-Tval:Tval/K:Tval-Tval/K)'; % Spatial mesh
        
        scfac1 = exp(-1i*k0*Tval/ep);
        scfac2 = dt*pi/Tval;
        cgsc = pi*cg/(ep*Tval);
        s = sign(k0);
        
        % Bright/Focusing initial conditions

        tfac = exp(1i*ad*dt/2);
        tfac2 = exp(1i*ad*dt);
        uint = sqrt(2*ad/anl)*sech(X);
     
        [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
        n2 = a1*uint.^2;
        w = fft(uint); 
        n0 = om*(2*Om*s-om)/(1+cg*om)*abs(uint).^2;        
        sint = ep^2*n0 + 2*ep*real(uint.*exp(1i*k0*X/ep) + ep*n2.*exp(2*1i*k0*X/ep));           
        
        Nvorts = 2;
        Lind = K;
        Rind = K+2;
        
        xint1 = pi/Tval*(X(Lind) + Tval);
        xint2 = pi/Tval*(X(Rind) + Tval);
        
        xpos = [xint1; xint2];
        zpos = [sint(Lind);sint(Rind)];
        
        xtrack = zeros(Nvorts,nmax+1);
        ztrack = zeros(Nvorts,nmax+1);
        
        sdrift = zeros(nmax+1,2);
                 
        sdrift(1,1) = Stokes_Drift_mean(Tval*xpos(1)/pi,Tval,w,K,s,Om,om,k0);
        sdrift(1,2) = Stokes_Drift_mean(Tval*xpos(2)/pi,Tval,w,K,s,Om,om,k0);
                
        xtrack(:,1) = (-Tval + Tval/pi*xpos)/ep;
        ztrack(:,1) = zpos;   
        
    % Solve NLS equation in time
    
        for nn=1:nmax
            
            t = (nn-1)*dt;
            af = w*tfac;
            cf = w*tfac2;
            
            for jj = 1:Nvorts
            
                theta = scfac1*exp(1i*(Tval*k0*xpos(jj)/pi + Om*t/ep)/ep);
                wc = ifft(w);
                xlocpos = Tval/pi*(xpos(jj)+cgsc*t) - Tval;
                xwave = Tval/pi*xpos(jj) - Tval;
                indc = loc_find(xlocpos,X);
                wval = wc(indc);
                n0 = om*(2*Om*s-om)/(1+cg*om)*abs(wval).^2;        
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*t/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*t/ep^2)) );
                
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*t)/pi,zval,Tval,w,K,s,cg,Om,om,k0,sig,ep);
                ax = scfac2*(om*zval/ep + xdot);
                
                theta = scfac1*exp(1i*(Tval*k0*(xpos(jj) + ax/2)/pi + Om*(t+dt/2)/ep)/ep);
                wph = ifft(af);
                xlocpos = Tval/pi*(xpos(jj)+cgsc*(t+dt/2)+ax/2) - Tval;
                xwave = Tval/pi*(xpos(jj)+ax/2) - Tval;
                indc = loc_find(xlocpos,X);
                wval = wph(indc);
                n0 = om*(2*Om*s-om)/(1+cg*om)*abs(wval).^2;        
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*(t+dt/2)/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*(t+dt/2)/ep^2)) );
                            
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*(t+dt/2)+ax/2)/pi,zval,Tval,af,K,s,cg,Om,om,k0,sig,ep);
                bx = scfac2*(om*zval/ep + xdot);
                
                theta = scfac1*exp(1i*(Tval*k0*(xpos(jj) + bx/2)/pi + Om*(t+dt/2)/ep)/ep);
                xlocpos = Tval/pi*(xpos(jj)+cgsc*(t+dt/2)+bx/2) - Tval;
                xwave = Tval/pi*(xpos(jj) + bx/2) - Tval;
                indc = loc_find(xlocpos,X);
                wval = wph(indc);
                n0 = om*(2*Om*s-om)/(1+cg*om)*abs(wval).^2;        
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*(t+dt/2)/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*(t+dt/2)/ep^2)) );
                
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*(t+dt/2)+bx/2)/pi,zval,Tval,af,K,s,cg,Om,om,k0,sig,ep);
                cx = scfac2*(om*zval/ep + xdot);
                
                theta = scfac1*exp(1i*(Tval*k0*(xpos(jj) + cx)/pi + Om*(t+dt)/ep)/ep);
                wpf = ifft(cf);
                xlocpos = Tval/pi*(xpos(jj)+cgsc*(t+dt)+cx) - Tval;
                xwave = Tval/pi*(xpos(jj) + cx) - Tval;
                indc = loc_find(xlocpos,X);
                wval = wpf(indc);
                n0 = om*(2*Om*s-om)/(1+cg*om)*abs(wval).^2;        
                n2 = a1*wval.^2;
                zval = ep^2*n0 + 2*ep*real(wval.*exp(1i*(k0*xwave/ep+ Om*(t+dt)/ep^2)) + ep*n2.*exp(2*1i*(k0*xwave/ep+ Om*(t+dt)/ep^2)) );
                
                xdot = phi_eval_1d(theta,Tval*(xpos(jj)+cgsc*(t+dt)+cx)/pi,zval,Tval,cf,K,s,cg,Om,om,k0,sig,ep);
                dx = scfac2*(om*zval/ep + xdot);
                
                xpos(jj) = xpos(jj) + (ax + 2*(bx+cx) + dx)/6;
                zpos(jj) = zval;
                
            end    
            
            w = w*tfac2;                
            
            xtrack(:,nn+1) = (-Tval + Tval*xpos/pi)/ep;
            ztrack(:,nn+1) = zpos;
            
            for jj=1:Nvorts
               sdrift(nn,jj) = Stokes_Drift_mean(Tval*(xpos(jj)+cgsc*(t+dt))/pi,Tval,w,K,s,Om,om,k0); 
            end
                     
        end
        
        w = ifft(w);
        no_cells = 1+mod(round(cg*(t+dt)/ep*K/Tval),2*K);
        w = w([no_cells:2*K 1:no_cells-1]');
        n0 = om*(2*Om*s-om)/(1+cg*om)*abs(w).^2;        
            
        n2 = a1*w.^2;
        w = ep^2*n0 + 2*ep*real(w.*exp(1i*(k0*X/ep + Om*(t+dt)/ep^2)) + ep*n2.*exp(2*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) );        
        
end