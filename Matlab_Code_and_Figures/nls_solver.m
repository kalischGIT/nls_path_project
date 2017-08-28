function [xtrack,ztrack,zdiff,w,sdrift] = nls_solver(K,Llx,ad,anl,cg,k0,Om,om,sig,ep,tf,dt)

    % solves u_t = dk2Om*u_xx + a3*|u|^2*u 

    % Llx is domain size, i.e. solve on domain [-Llx,Llx]
    % sig is -1 for bright/focusing, +1 for dark/defocusing
    % tf is total simulation time
    % mu is magnitude of damping on envelope of dark initial condition
   
        nmax = round(tf/dt); % Step size for time integrator and number of time steps
    
        Dx = 1i*pi/Llx*[(0:K-1) 0 (-K+1:-1)].'; % Fourier transform of first derivative
        
        Dx2 = Dx.^2; % Fourier transform of second derivative
            
        X = (-Llx:Llx/K:Llx-Llx/K)'; % Spatial mesh
        
        scfac1 = exp(-1i*k0*Llx/ep);
        scfac2 = dt*pi/Llx;
        cgsc = pi*cg/(ep*Llx);
        
        % Bright/Focusing initial conditions
    
        uint = sqrt(2*real(ad/anl))*sech(X); 
        [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
        n2 = a1*uint.^2;
        n3 = a3*uint.^3;
        
        sint = 2*ep*real(uint.*exp(1i*k0*X/ep) + ep*n2.*exp(2*1i*k0*X/ep) + ep^2*n3.*exp(3*1i*k0*X/ep));
                
        w = fft(uint); 
            
        Nvorts = 2;
        Lind = K/4;
        Rind = 3*K/4;
        
        xint1 = pi/Llx*(X(Lind) + Llx);
        xint2 = pi/Llx*(X(Rind) + Llx);
        
        xpos = [xint1; xint2];
        zpos = [sint(Lind);sint(Rind)];
        
        xtrack = zeros(Nvorts,nmax+1);
        ztrack = zeros(Nvorts,nmax+1);
        zdiff = zeros(Nvorts,nmax+1);
        sdrift = zeros(nmax+1,2);
        
        fuint = fft(uint);
        
        sdrift(1,1) = log10(abs(Stokes_Drift(Llx*xpos(1)/pi,zpos(1),Llx,fuint,K,sign(k0),Om,k0,ep)));
        sdrift(1,2) = log10(abs(Stokes_Drift(Llx*xpos(2)/pi,zpos(2),Llx,fuint,K,sign(k0),Om,k0,ep)));
        
        xtrack(:,1) = (-Llx + Llx/pi*xpos)/ep;
        ztrack(:,1) = zpos;       
        
    % Begin setup of RK4 method

        Lap = dt*ad*Dx2/2;
        E = exp(Lap);
        E2 = E.^2;
                
    % Solve NLS equation in time
    
        for nn=1:nmax
            
            t = (nn-1)*dt;
            
            wp = ifft(w);
            a = dt*anl*fft(wp.^2.*conj(wp));
                        
            af = E.*(w+a/2);
            ap = ifft(af);
            b = dt*anl*fft(ap.^2.*conj(ap));
            
            bf = E.*w + b/2;
            bp = ifft(bf);
            c = dt*anl*fft(bp.^2.*conj(bp));
            
            cf = E2.*w + E.*c;
            cp = ifft(cf);
            d = dt*anl*fft(cp.^2.*conj(cp));
            
            for jj = 1:Nvorts
            
                theta = scfac1*exp(1i*(Llx*k0*xpos(jj)/pi + Om*t/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Llx*(xpos(jj)+cgsc*t)/pi,zpos(jj),Llx,w,K,sign(k0),cg,Om,om,k0,sig,ep);
                ax = scfac2*(-om*zpos(jj)/ep + xdot);
                az = dt*zdot/ep;            
                        
                theta = scfac1*exp(1i*(Llx*k0*(xpos(jj) + ax/2)/pi + Om*(t+dt/2)/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Llx*(xpos(jj)+cgsc*(t+dt/2)+ax/2)/pi,zpos(jj)+az/2,Llx,af,K,sign(k0),cg,Om,om,k0,sig,ep);
                bx = scfac2*(-om*(zpos(jj)+az/2)/ep + xdot);
                bz = dt*zdot/ep;            
                       
                theta = scfac1*exp(1i*(Llx*k0*(xpos(jj) + bx/2)/pi + Om*(t+dt/2)/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Llx*(xpos(jj)+cgsc*(t+dt/2)+bx/2)/pi,zpos(jj)+bz/2,Llx,bf,K,sign(k0),cg,Om,om,k0,sig,ep);
                cx = scfac2*(-om*(zpos(jj)+bz/2)/ep + xdot);
                cz = dt*zdot/ep;            
                        
                theta = scfac1*exp(1i*(Llx*k0*(xpos(jj) + cx)/pi + Om*(t+dt)/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Llx*(xpos(jj)+cgsc*(t+dt)+cx)/pi,zpos(jj)+cz,Llx,cf,K,sign(k0),cg,Om,om,k0,sig,ep);
                dx = scfac2*(-om*(zpos(jj)+cz)/ep + xdot);
                dz = dt*zdot/ep;            
                
                xpos(jj) = xpos(jj) + (ax + 2*(bx+cx) + dx)/6;
                zpos(jj) = zpos(jj) + (az + 2*(bz+cz) + dz)/6;                  
                
            end    
            
            w = E2.*w + (E2.*a + 2*E.*(b+c) + d)/6;                
            
            for jj=1:Nvorts
               sdrift(nn,jj) = log10(abs(Stokes_Drift(Llx*(xpos(jj)+cgsc*(t+dt))/pi,zpos(jj),Llx,w,K,sign(k0),Om,k0,ep))); 
            end
            
            xtrack(:,nn+1) = (-Llx + Llx*xpos/pi)/ep;
            ztrack(:,nn+1) = zpos;
            
            wp = ifft(w);
            
            no_cells = 1+mod(round(cg*(t+dt)/ep*K/Llx),2*K);
            wp = wp([no_cells:2*K 1:no_cells-1]');
            
            n2 = a1*wp.^2;
            n3 = a3*wp.^2;
            wp = 2*ep*real( wp.*exp(1i*(k0*X/ep+ Om*(t+dt)/ep^2)) + ep*n2.*exp(2*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) + ep^2*n3.*exp(3*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) );
            
            for ll=1:Nvorts
                yy = spline(X/ep,wp,xtrack(ll,nn+1));
                zdiff(ll,nn+1) = log10(abs(zpos(ll)-yy));
            end
            
        end
        
        w = ifft(w);
        no_cells = 1+mod(round(cg*t/ep*K/Llx),2*K);
        w = w([no_cells:2*K 1:no_cells-1]');
        
        n2 = a1*w.^2;
        n3 = a3*w.^3;
        w = 2*ep*real( w.*exp(1i*(k0*X/ep + Om*(t+dt)/ep^2)) + ep*n2.*exp(2*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) + ep^2*n3.*exp(3*1i*(k0*X/ep+ Om*(t+dt)/ep^2)));        
        
end