function nls_path_plotter_periodic_1d(ep,kval,K,tf)

    sig = 1e-5;
    Tval = 2*ellipke(kval^2);
    X = (-Tval:Tval/K:Tval-Tval/K)'; % Spatial mesh
    
    mval = round(Tval/(pi*ep));
    k0 = pi*mval*ep/Tval;
    
    dt = 5e-4;
    nmax = round(tf/dt);
    tvals = 1/ep^2*linspace(0,tf,nmax+1);
    %omval = zero_lan(k0,1e-3,1e-13);
    %disp(omval)
    omval = 1.12;
    
    omvals = linspace(-20,20,2e3)';
    ampvals = zeros(length(omvals),1);
    stksvals = zeros(length(omvals),1);

    for jj=1:length(omvals)
        [Om,~,ad,anl] = param_maker(k0,omvals(jj),sig);
        ampvals(jj) = sqrt(2*abs(ad/anl));
        stksvals(jj) = -2*k0*Om*(2-sign(k0)*omvals(jj)/Om)*ampvals(jj)^2;
    end
    
    %[val,ind] = max(ampvals);
        
    % Zero case
    [Om,cg,ad,anl] = param_maker(k0,0,sig);
    [xtrack0,ztrack0,xmtrack0,zmtrack0,u0,sdrift0] = nls_stepper_periodic_1d(K,ad,anl,cg,k0,kval,Om,0,sig,ep,tf,dt);
    disp(sign(ad/anl))
    [Om,cg,ad,anl] = param_maker(k0,omval,sig);
    [xtrack1,ztrack1,xmtrack1,zmtrack1,u4,sdrift1] = nls_stepper_periodic_1d(K,ad,anl,cg,k0,kval,Om,omval,sig,ep,tf,dt);
    disp(sign(ad/anl))
    [Om,cg,ad,anl] = param_maker(k0,-1,sig);
    [xtrackn1,ztrackn1,xmtrackn1,zmtrackn1,~,sdriftn1] = nls_stepper_periodic_1d(K,ad,anl,cg,k0,kval,Om,-1,sig,ep,tf,dt);
    disp(sign(ad/anl))
    clf
    
    figure(1)
    hold on
        
    plot(xtrack1(2,1:end),ztrack1(2,1:end),'k','LineWidth',2)
    plot(xtrack1(2,1),ztrack1(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    plot(xtrack1(2,end),ztrack1(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
    
    plot(X/ep,u4,'k','LineWidth',2)
    
    hold off
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    
    figure(2)
    
    hold on
        
    %plot(xmtrack0(2,1:end),zmtrack0(2,1:end),'k','LineWidth',2)
    %plot(xmtrack0(2,1),zmtrack0(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    %plot(xmtrack0(2,end),zmtrack0(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
    
    plot(xtrack0(2,1:end),ztrack0(2,1:end),'k','LineWidth',2)
    plot(xtrack0(2,1),ztrack0(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    plot(xtrack0(2,end),ztrack0(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
        
    hold off
        
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    
    figure(3)
    
    hold on
        
    %plot(xmtrack1(2,1:end),zmtrack1(2,1:end),'k','LineWidth',2)
    %plot(xmtrack1(2,1),zmtrack1(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    %plot(xmtrack1(2,end),zmtrack1(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
    
    plot(xtrack1(2,1:end),ztrack1(2,1:end),'k','LineWidth',2)
    plot(xtrack1(2,1),ztrack1(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    plot(xtrack1(2,end),ztrack1(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
    
    hold off
        
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    
    figure(4)
    
    hold on
        
    %plot(xmtrackn1(2,1:end),zmtrackn1(2,1:end),'k','LineWidth',2)
    %plot(xmtrackn1(2,1),zmtrackn1(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    %plot(xmtrackn1(2,end),zmtrackn1(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
    
    plot(xtrackn1(2,1:end),ztrackn1(2,1:end),'k','LineWidth',2)
    plot(xtrackn1(2,1),ztrackn1(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    plot(xtrackn1(2,end),ztrackn1(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
        
    hold off
        
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    
    figure(5)
    hold on 
    plot(tvals(1:end-2),sdriftn1(1:end-2,2),'k--','LineWidth',2)
    plot(tvals(1:end-2),sdrift1(1:end-2,2),'k-.','LineWidth',2)
    plot(tvals(1:end-2),sdrift0(1:end-2,2),'k','LineWidth',2)    
    hold off
    
    haxes = set(gca,'FontSize',30);
    set(haxes,'Interpreter','LaTeX')
    hlegend = legend('$\omega=-1$','$\omega=1$','$\omega=0$');
    set(hlegend,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX')
    ylabel('$\tilde{u}^{S}(x_{m}(t),0,t)$','Interpreter','LaTeX')  
    
    figure(6)
    plot(omvals,ampvals,'k','LineWidth',2)
    haxes = set(gca,'FontSize',30);
    set(haxes,'Interpreter','LaTeX')
    xlabel('$\omega$','Interpreter','LaTeX')
    ylabel('$\sqrt{2 \frac{\alpha_{d}}{\alpha_{nl}}}$','Interpreter','LaTeX')        
    
    figure(7)
    plot(omvals,stksvals,'k','LineWidth',2)
    haxes = set(gca,'FontSize',30);
    set(haxes,'Interpreter','LaTeX')
    xlabel('$\omega$','Interpreter','LaTeX')
    ylabel('$\tilde{u}^{S}$','Interpreter','LaTeX')           