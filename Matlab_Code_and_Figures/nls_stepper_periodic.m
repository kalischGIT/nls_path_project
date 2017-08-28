function [xtrack,ztrack,xmtrack,zmtrack,zdiff,w,sdrift] = nls_stepper_periodic(K,ad,anl,cg,k0,kval,Om,om,sig,ep,tf,dt)

    % Llx is domain size, i.e. solve on domain [-Llx,Llx]
    % sig is -1 for bright/focusing, +1 for dark/defocusing
    % tf is total simulation time
    % mu is magnitude of damping on envelope of dark initial condition
        nmax = round(tf/dt); % Step size for time integrator and number of time steps
        Tval = 2*ellipke(kval^2);
        X = (-Tval:Tval/K:Tval-Tval/K)'; % Spatial mesh
        Dx = pi*1i/Tval*[0:K -K+1:-1]';
        
        scfac1 = exp(-1i*k0*Tval/ep);
        scfac2 = dt*pi/Tval;
        cgsc = pi*cg/(ep*Tval);
        s = sign(k0);
        
        % Bright/Focusing initial conditions

        if ad/anl > 0
            [~,cnvals,~] = ellipj(X,kval^2);
            uint = kval*sqrt(2*ad/anl)*cnvals; 
            tfac = exp(-1i*ad*(1-2*kval^2)*dt/2);
            tfac2 = exp(-1i*ad*(1-2*kval^2)*dt);            
        else
            [snvals,~,~] = ellipj(X,kval^2);
            uint = kval*sqrt(-2*ad/anl)*snvals; 
            tfac = exp(-1i*ad*(1+kval^2)*dt/2);
            tfac2 = exp(-1i*ad*(1+kval^2)*dt);
        end
        %{
            f = @(x) kval*sqrt(2*ad/anl)*(sech(x) - (1-kval^2)*(sinh(x).*cosh(x)-x).*tanh(x).*sech(x));
            plot(X,uint,'k')
            pause
        %}
        [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig);
    
        n2 = a1*uint.^2;
        n3 = a3*uint.^3;
        w = fft(uint); 
        q0 = -(om-2*s*Om)*real(ifft(abs(Dx).*fft(abs(uint).^2)));
        n0 = -cg/(1+cg*om)*q0;        
        sint = ep^3*n0 + 2*ep*real(uint.*exp(1i*k0*X/ep) + ep*n2.*exp(2*1i*k0*X/ep) + ep^2*n3.*exp(3*1i*k0*X/ep));           
        
        Nvorts = 2;
        Lind = K;
        Rind = K+2;
        
        xint1 = pi/Tval*(X(Lind) + Tval);
        xint2 = pi/Tval*(X(Rind) + Tval);
        
        xpos = [xint1; xint2];
        zpos = [sint(Lind);sint(Rind)];
        
        xmpos = [xint1; xint2];
        zmpos = [0;0];
        %y20 = [sint(Lind)-zmpos(1);sint(Rind)-zmpos(2)];
        
        
        xtrack = zeros(Nvorts,nmax+1);
        ztrack = zeros(Nvorts,nmax+1);
        
        xmtrack = zeros(Nvorts,nmax+1);
        zmtrack = zeros(Nvorts,nmax+1);
        
        zdiff = zeros(Nvorts,nmax+1);
        sdrift = zeros(nmax+1,2);
                 
        %sdrift(1,1) = Stokes_Drift(Tval*xpos(1)/pi,Tval,w,K,s,Om,om,k0);
        %sdrift(1,2) = Stokes_Drift(Tval*xpos(2)/pi,Tval,w,K,s,Om,om,k0);
        
        sdrift(1,1) = Stokes_Drift_mean(Tval*xmpos(1)/pi,Tval,w,K,s,Om,om,k0);
        sdrift(1,2) = Stokes_Drift_mean(Tval*xmpos(2)/pi,Tval,w,K,s,Om,om,k0);
                
        xtrack(:,1) = (-Tval + Tval/pi*xpos)/ep;
        ztrack(:,1) = zpos;   
        
        xmtrack(:,1) = (-Tval + Tval/pi*xmpos)/ep;
        zmtrack(:,1) = zmpos;   
                
    % Solve NLS equation in time
    
        for nn=1:nmax
            
            t = (nn-1)*dt;
            af = w*tfac;
            cf = w*tfac2;
            
            for jj = 1:Nvorts
            
                theta = scfac1*exp(1i*(Tval*k0*xpos(jj)/pi + Om*t/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Tval*(xpos(jj)+cgsc*t)/pi,zpos(jj),Tval,w,K,s,cg,Om,om,k0,sig,ep);
                [xmdot,zmdot] = phi_eval_mean(Tval*(xmpos(jj)+cgsc*t)/pi,zmpos(jj),Tval,w,K,s,cg,Om,om,k0,sig,ep);
                
                ax = scfac2*(om*zpos(jj)/ep + xdot);
                az = dt*zdot/ep;            
                
                amx = scfac2*(om*zmpos(jj)/ep + ep*xmdot);
                amz = dt*zmdot;
                
                theta = scfac1*exp(1i*(Tval*k0*(xpos(jj) + ax/2)/pi + Om*(t+dt/2)/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Tval*(xpos(jj)+cgsc*(t+dt/2)+ax/2)/pi,zpos(jj)+az/2,Tval,af,K,s,cg,Om,om,k0,sig,ep);
                [xmdot,zmdot] = phi_eval_mean(Tval*(xmpos(jj)+cgsc*(t+dt/2)+amx/2)/pi,zmpos(jj)+amz/2,Tval,w,K,s,cg,Om,om,k0,sig,ep);
                
                bx = scfac2*(om*(zpos(jj)+az/2)/ep + xdot);
                bz = dt*zdot/ep;            
                
                bmx = scfac2*(om*(zmpos(jj)+amz/2)/ep + ep*xmdot);
                bmz = dt*zmdot;
                       
                theta = scfac1*exp(1i*(Tval*k0*(xpos(jj) + bx/2)/pi + Om*(t+dt/2)/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Tval*(xpos(jj)+cgsc*(t+dt/2)+bx/2)/pi,zpos(jj)+bz/2,Tval,af,K,s,cg,Om,om,k0,sig,ep);
                [xmdot,zmdot] = phi_eval_mean(Tval*(xmpos(jj)+cgsc*(t+dt/2)+bmx/2)/pi,zmpos(jj)+bmz/2,Tval,w,K,s,cg,Om,om,k0,sig,ep);
                                
                cx = scfac2*(om*(zpos(jj)+bz/2)/ep + xdot);
                cz = dt*zdot/ep;            
                
                cmx = scfac2*(om*(zmpos(jj)+bmz/2)/ep + ep*xmdot);
                cmz = dt*zmdot;
                        
                theta = scfac1*exp(1i*(Tval*k0*(xpos(jj) + cx)/pi + Om*(t+dt)/ep)/ep);
                [xdot,zdot] = phi_eval(theta,Tval*(xpos(jj)+cgsc*(t+dt)+cx)/pi,zpos(jj)+cz,Tval,cf,K,s,cg,Om,om,k0,sig,ep);
                [xmdot,zmdot] = phi_eval_mean(Tval*(xmpos(jj)+cgsc*(t+dt)+cmx)/pi,zmpos(jj)+cmz,Tval,w,K,s,cg,Om,om,k0,sig,ep);
                
                dx = scfac2*(om*(zpos(jj)+cz)/ep + xdot);
                dz = dt*zdot/ep;            
                
                dmx = scfac2*(om*(zmpos(jj)+cmz)/ep + ep*xmdot);
                dmz = dt*zmdot;
                
                xpos(jj) = xpos(jj) + (ax + 2*(bx+cx) + dx)/6;
                zpos(jj) = zpos(jj) + (az + 2*(bz+cz) + dz)/6;                  
                
                xmpos(jj) = xmpos(jj) + (amx + 2*(bmx+cmx) + dmx)/6;
                zmpos(jj) = zmpos(jj) + (amz + 2*(bmz+cmz) + dmz)/6;                  
                 
            end    
            
            w = w*tfac2;                
            
            xtrack(:,nn+1) = (-Tval + Tval*xpos/pi)/ep;
            ztrack(:,nn+1) = zpos;
            
            xmtrack(:,nn+1) = (-Tval + Tval*xmpos/pi)/ep;
            zmtrack(:,nn+1) = zmpos;
                        
            wp = ifft(w);
            
            no_cells = 1+mod(round(cg*(t+dt)/ep*K/Tval),2*K);
            wp = wp([no_cells:2*K 1:no_cells-1]');
            %n0 = cg*(om-2*s*Om)/(1+cg*(1+om))*real(ifft(abs(Dx).*fft(abs(wp).^2)));        
            n0 = cg*(om-2*s*Om)/(1+cg*om)*real(ifft(abs(Dx).*fft(abs(wp).^2)));        
            n2 = a1*wp.^2;
            n3 = a3*wp.^2;
            wp = ep^3*n0 + 2*ep*real(wp.*exp(1i*(k0*X/ep+ Om*(t+dt)/ep^2)) + ep*n2.*exp(2*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) + ep^2*n3.*exp(3*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) );
            
            for jj=1:Nvorts
               %sdrift(nn,jj) = Stokes_Drift(Tval*(xpos(jj)+cgsc*(t+dt))/pi,Tval,w,K,s,Om,om,k0); 
               sdrift(nn,jj) = Stokes_Drift_mean(Tval*(xmpos(jj)+cgsc*(t+dt))/pi,Tval,w,K,s,Om,om,k0); 
            end
            
            for ll=1:Nvorts
                yy = spline(X/ep,wp,xtrack(ll,nn+1));
                zdiff(ll,nn+1) = (abs(zpos(ll)-yy));
            end
            
        end
        
        w = ifft(w);
        no_cells = 1+mod(round(cg*t/ep*K/Tval),2*K);
        w = w([no_cells:2*K 1:no_cells-1]');
        %n0 = cg*(om-2*s*Om)/(1+cg*(1+om))*real(ifft(abs(Dx).*fft(abs(w).^2)));        
        n0 = cg*(om-2*s*Om)/(1+cg*om)*real(ifft(abs(Dx).*fft(abs(w).^2)));        
        
        n2 = a1*w.^2;
        n3 = a3*w.^3;
        w = ep^3*n0 + 2*ep*real(w.*exp(1i*(k0*X/ep + Om*(t+dt)/ep^2)) + ep*n2.*exp(2*1i*(k0*X/ep+ Om*(t+dt)/ep^2)) + ep^2*n3.*exp(3*1i*(k0*X/ep+ Om*(t+dt)/ep^2)));        
        
end