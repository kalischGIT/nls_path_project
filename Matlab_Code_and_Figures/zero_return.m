function zero_return(Pvals,sig,tol)

kvals = linspace(0,20,Pvals);
omvals = zeros(length(kvals),1);
%omvals_wrong = zeros(length(kvals),1);

for jj=1:length(kvals)
    omvals(jj) = zero_lan(kvals(jj),sig,tol);
    %omvals_wrong(jj) = zero_lan_wrong(kvals(jj),sig,tol);
end

figure(1)
%plot(kvals,omvals,'k-',kvals,omvals_wrong,'k--','LineWidth',2)
plot(kvals,omvals,'k-','LineWidth',2)
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
xlabel('$k_{0}$','Interpreter','LaTeX','FontSize',30)
ylabel('$\omega$','Interpreter','LaTeX','FontSize',30)
%legend({'$u^{L}_{p}(k_{0},\omega)$','$u^{L,a}_{p}(k_{0},\omega)$'},'Interpreter','LaTeX')