function g = gfun_nl_oscillator(mp,ms,kp,ks,zetap,zetas,Fs,S0)

omegap = (kp./mp).^0.5;
omegas = (ks./ms).^0.5;

gamma  = ms./mp;

omegaa = (omegap+omegas)/2;
zetaa  = (zetap+zetas)/2;

theta  = (omegap-omegas)./omegaa;

Exs2   = pi*S0./(4*zetas.*omegas.^3) .* ( ...
   zetaa.*zetas ./ ( zetap.*zetas.*(4*zetaa.^2+theta.^2) + gamma.*zetaa.^2 ) .* ...
   ( (zetap.*omegap.^3+zetas.*omegas.^3).*omegap ./ (4*zetaa.*omegaa.^4) ) ...
                                  );
g = Fs - ks*3.*Exs2.^0.5;