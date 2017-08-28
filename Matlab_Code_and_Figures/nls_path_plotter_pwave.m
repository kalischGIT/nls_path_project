function nls_path_plotter_pwave(ep,K,Llx,tf)

    sig = 1e-3;
    
    k0 = 1;
    %omval = 1.682810828831611;
    omval = 4;
    A = 1;
    
    omvals = linspace(-omval,omval,1e3);
    anlvals = zeros(length(omvals),1);
    advals = zeros(length(omvals),1);
    cgvals = zeros(length(omvals),1);
    
    for jj=1:length(omvals)
        [~,cgvals(jj),advals(jj),anlvals(jj)] = param_maker(k0,omvals(jj),sig);
    end
    
    [Om,cg,ad,anl] = param_maker(k0,omval,sig);
    [Tvals,xtrack0,ztrack0,u0,sdrift0] = nls_solver_pwave_1d(K,Llx,A,cg,anl,k0,Om,omval,sig,ep,tf);
    disp(sign(ad/anl))  
    disp(min(sdrift0))
    disp(max(sdrift0))
    clf

    figure(1)
    hold on
        
    plot(xtrack0,ztrack0,'k','LineWidth',2)
    plot(xtrack0(1),ztrack0(1),'.','MarkerSize',58','color',[0.8 0.8 0.8]);
    plot(xtrack0(end),ztrack0(end),'.','MarkerSize',58','color',[0.1 0.1 0.1]);
        
    hold off
        
    haxes = set(gca,'FontSize',30);
    set(haxes,'Interpreter','LaTeX')
    xlabel('$x$','Interpreter','LaTeX')    
    ylabel('$z$','Interpreter','LaTeX')  
    
    figure(2)
    plot(omvals,cgvals,'k','LineWidth',2)
    haxes = set(gca,'FontSize',30);
    set(haxes,'Interpreter','LaTeX')
    xlabel('$\omega$','Interpreter','LaTeX')    
    ylabel('$c_{g}(1,\omega)$','Interpreter','LaTeX')  