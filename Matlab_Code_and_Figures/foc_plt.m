function foc_plt(sig)

Om = @(k0,om) .5*(sign(k0).*om + sqrt(om.^2 + 4.*abs(k0).*(1+sig.*k0.^2)));

cg = @(k0,om) (1+3*sig*k0.^2)./(2*sign(k0).*Om(k0,om)-om);

ad = @(k0,om) ((cg(k0,om)).^2 - 3*abs(k0)*sig)./(2*Om(k0,om) - sign(k0).*om);

av = @(k0,om) (cg(k0,om).*abs(k0) - 2*Om(k0,om).*sign(k0)).*om.^4 + k0.*(4*k0.^2.*sign(k0).*sig + 2*Om(k0,om).*cg(k0,om) - sign(k0)).*om.^3 ...
     + k0.*(16*cg(k0,om).*k0.^3*sig - 8*Om(k0,om).*k0.^2*sig + 10*cg(k0,om).*k0 - 6*Om(k0,om)).*om.^2 ...
     - k0.^2.*(15*Om(k0,om).*cg(k0,om).*k0.^2.*sign(k0)*sig - 16*k0.^4*sig^2 - 24*k0.^2*sig - 2).*om ...
     + k0.^3.*(2*cg(k0,om)*k0.^4.*sign(k0)*sig^2 + cg(k0,om).*k0.^2.*sign(k0)*sig - 15*Om(k0,om).*k0.*sign(k0)*sig + 8*cg(k0,om).*sign(k0));

   
anl = @(k0,om) k0.*(sign(k0).*(k0.^3).*(8 + sig*k0.^2 + 2*(sig*k0.^2).^2)+om.*av(k0,om))./...
              ((2.*sign(k0).*Om(k0,om)-om).*(1+cg(k0,om).*om).*(4*(Om(k0,om)).^2-sign(k0).*(2.*k0.*(1+4.*sig.*k0.^2)+2.*om.*Om(k0,om))));

n0 = @(k0,om) om.^2.*(2.*Om(k0,om).*sign(k0)-om)./(1+om.*cg(k0,om));
us = @(k0,om) -2.*k0.*Om(k0,om).*(2 - 2*sign(k0).*om./Om(k0,om) + (om./Om(k0,om)).^2);          
          
afoc = @(k0,om) anl(k0,om)./ad(k0,om);

kvals = linspace(-20,20,2e3);
ovals = linspace(-40,40,2e3);
tvals = zeros(length(ovals),length(kvals));
ampvals = zeros(length(ovals),length(kvals));
stksvals = zeros(length(ovals),length(kvals));
lansvals = zeros(length(ovals),length(kvals));

for jj=1:length(ovals)
    tvals(:,jj) = sign(afoc(kvals(jj),ovals));
    ampvals(:,jj) = log10(sqrt(2*abs(ad(kvals(jj),ovals)./anl(kvals(jj),ovals))));    
    stksvals(:,jj) = log10( abs( 4*abs(ad(kvals(jj),ovals)./anl(kvals(jj),ovals)).*kvals(jj).*Om(kvals(jj),ovals).*(2 - 2*sign(kvals(jj)).*ovals./Om(kvals(jj),ovals) + (ovals./Om(kvals(jj),ovals)).^2) )) ;
    lansvals(:,jj) = log10(abs( ( n0(kvals(jj),ovals) + us(kvals(jj),ovals) ) ).*2.*abs( ad(kvals(jj),ovals)./anl(kvals(jj),ovals) ) );
end

figure(1)
surf(kvals,ovals,tvals,'LineStyle','none')
view(2)
haxes = set(gca,'FontSize',30);
set(haxes,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX','FontSize',30)    
ylabel('$\omega$','Interpreter','LaTeX','FontSize',30)  

figure(2)
surf(kvals,ovals,ampvals,'LineStyle','none')
view(2)
haxes = set(gca,'FontSize',30);
set(haxes,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX','FontSize',30)    
ylabel('$\omega$','Interpreter','LaTeX','FontSize',30)  

figure(3)
surf(kvals,ovals,stksvals,'LineStyle','none')
view(2)
haxes = set(gca,'FontSize',30);
set(haxes,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX','FontSize',30)    
ylabel('$\omega$','Interpreter','LaTeX','FontSize',30)  

figure(4)
surf(kvals,ovals,lansvals,'LineStyle','none')
view(2)
haxes = set(gca,'FontSize',30);
set(haxes,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX','FontSize',30)    
ylabel('$\omega$','Interpreter','LaTeX','FontSize',30)  

%{
figure(2)
[C,h] = contour(kvals,ovals,log10(ampvals),200,'k');
vlevels = [-3,0,2];
clabel(C,h,vlevels,'FontSize',30,'Color','black','LabelSpacing',100)
haxes = set(gca,'FontSize',30);
set(haxes,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX')    
ylabel('$\omega$','Interpreter','LaTeX')  
%}