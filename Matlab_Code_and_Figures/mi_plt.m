function mi_plt(sig)

Omp = @(k0,om) .5*(-sign(k0).*om + sqrt(om.^2 + 4.*abs(k0).*(1+sig.*k0.^2)));
Omn = @(k0,om) .5*(-sign(k0).*om - sqrt(om.^2 + 4.*abs(k0).*(1+sig.*k0.^2)));

kvals = linspace(-20,20,4e3);
ovals = linspace(-40,40,4e3);
pvals = zeros(length(ovals),length(kvals));
nvals = zeros(length(ovals),length(kvals));

for jj=1:length(kvals)
    pvals(:,jj) = sign(Omp(kvals(jj),ovals));    
    nvals(:,jj) = sign(Omn(kvals(jj),ovals));    
end

figure(1)
surf(kvals,ovals,pvals,'LineStyle','none')
view(2)
haxes = set(gca,'FontSize',30);
set(haxes,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX')    
ylabel('$\omega$','Interpreter','LaTeX')  

figure(2)
surf(kvals,ovals,nvals,'LineStyle','none')
view(2)
haxes = set(gca,'FontSize',30);
set(haxes,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX')    
ylabel('$\omega$','Interpreter','LaTeX')  