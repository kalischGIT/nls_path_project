function nls_path_plotter_one_at_a_time(ep,k0,kval,K,Llx,tf)

    sig = 1e-3;
    Tval = 2*Llx*ellipke(kval^2);
    X = (-Tval:Tval/K:Tval-Tval/K)'; % Spatial mesh
    disp(Tval)
    
    dt = 5e-4;
    nmax = round(tf/dt);
    tvals = 1/ep^2*linspace(0,tf,nmax+1);
    omval = 1;
    
    [Om,cg,ad,anl] = param_maker(k0,omval,sig);
    [xtrack1,ztrack1,zdiff1,u1,sdrift1] = nls_stepper(K,Llx,ad,anl,cg,k0,kval,Om,omval,sig,ep,tf,dt);
    disp(sign(ad/anl))
    disp(sign(Om))
    disp(anl)
    disp(sqrt(2*abs(ad/anl)))
    clf
    
    figure(1)
    hold on
        
    plot(xtrack1(2,1:end),ztrack1(2,1:end),'k','LineWidth',2)
    plot(xtrack1(2,1),ztrack1(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    plot(xtrack1(2,end),ztrack1(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
    
    plot(X/ep,u1,'k','LineWidth',2)
    
    hold off
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    
    figure(2)
    
    hold on
        
    plot(xtrack1(2,1:end),ztrack1(2,1:end),'k','LineWidth',2)
    plot(xtrack1(2,1),ztrack1(2,1),'.','MarkerSize',40','color',[0.8 0.8 0.8]);
    plot(xtrack1(2,end),ztrack1(2,end),'.','MarkerSize',40','color',[0.1 0.1 0.1]);
    
    hold off
        
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    
    figure(3)
    plot(tvals(1:end-2),sdrift1(1:end-2,2),'k-.','LineWidth',2)
    
    haxes = set(gca,'FontSize',30);
    set(haxes,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX')
    ylabel('$\tilde{u}^{S}(x_{s}(t),t)$','Interpreter','LaTeX')  
    
    figure(4)
    plot(tvals(2:end),zdiff1(2,2:end))