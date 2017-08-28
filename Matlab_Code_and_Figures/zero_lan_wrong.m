function w2 = zero_lan_wrong(k0,sig,tol)

s = sign(k0);

Om = @(w) .5*(s.*w + sqrt(w.^2 + 4*abs(k0)*(1+sig*k0^2)));
cg = @(w) (1+3*sig*k0^2)./(2*s*Om(w)-w);

f = @(w) w^2*(2*Om(w)*s-w)/(1+w*cg(w)) - 2*Om(w)*k0*(2 - s*w./Om(w));

w1 = .1;
w0 = -.05;

f1 = f(w1);
f0 = f(w0);

w2 = w1 - f1*(w1-w0)/(f1-f0);

while(abs(w2-w1)>=tol);
    
    w0 = w1;
    w1 = w2;
    f0 = f1;
    f1 = f(w1);

    w2 = w1 - f1*(w1-w0)/(f1-f0);

end

