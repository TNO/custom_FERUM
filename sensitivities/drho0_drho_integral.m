function drho0_drho = drho0_drho_integral(rho0,margi,margj,Z1,Z2,X1,X2,WIP,detJ)

% Compute the bi-folded drho/drho0-integral by 2D numerical integration, then takes its inverse

dPHI2_drho0 = 1/2/pi/(1-rho0^2)^(3/2)*exp(-1/(2-2*rho0^2)*(Z1.^2-2*rho0*Z1.*Z2+Z2.^2))*rho0 + 1/2/pi/(1-rho0^2)^(1/2)*(-4/(2-2*rho0^2)^2*(Z1.^2-2*rho0*Z1.*Z2+Z2.^2)*rho0+2/(2-2*rho0^2)*Z1.*Z2).*exp(-1/(2-2*rho0^2)*(Z1.^2-2*rho0*Z1.*Z2+Z2.^2));
drho_drho0 = sum ( sum ( (X1-margi(2))/margi(3) .* (X2-margj(2))/margj(3) .* dPHI2_drho0 * detJ .* WIP ) );

% Compute drho0/drho (inverse of drho/drho0)
drho0_drho = 1/drho_drho0;
