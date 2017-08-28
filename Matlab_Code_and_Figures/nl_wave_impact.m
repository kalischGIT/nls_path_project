function [omvals,ampvals] = nl_wave_impact(k0,sig)

omvals = linspace(-10,10,1e3)';
ampvals = zeros(length(omvals),1);

for jj=1:length(omvals)
    [Om,cg,ad,anl] = param_maker(k0,omvals(jj),sig);
    ampvals(jj) = sqrt(2*abs(ad/anl));
end