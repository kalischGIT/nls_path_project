function [a1,a2,a3,b3] = nls_expan_params(k0,Om,om,sig)
s = sign(k0);

a1 = (abs(k0)*(2*Om^2-2*s*om*Om+om^2) -2*k0*om*Om)/(4*Om^2-2*s*(k0*(1+4*sig*k0^2)+om*Om));
a2 = (k0*om*(2*k0*(1+4*sig*k0^2)+2*om*Om)-2*Om*k0*(2*Om^2-2*s*om*Om+om^2))/(4*Om^2-2*s*(k0*(1+4*sig*k0^2)+om*Om));

a3t = -1/2*(3*k0^5*sig - 2*Om^2*k0^2*s - 2*Om*a1*abs(k0)*om + 4*Om^2*a1*k0 - 2*Om*a2*abs(k0) + 2*Om*k0^2*om + 2*a1*k0*om^2 + 2*a2*k0*om);
b3t = -(2*Om*a1*abs(k0) - k0^2*s*om + a2*k0);
    
a3 = (s*a3t/3 - Om*b3t)/(s*(9*k0^3*sig+k0+om*Om)/3-Om^2);
b3 = (-Om*a3t+(9*k0^3*sig+k0+om*Om)*b3t)/(s*(9*k0^3*sig+k0+om*Om)/3-Om^2);
    