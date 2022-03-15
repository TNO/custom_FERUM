function f = abs_of_betacorr(rho0,rho_target,margi,margj,Z1,Z2,X1,X2,WIP,detJ)

% Compute the absolute value of the bi-folded rho-integral by 2D numerical integration

f = abs( rho_target - rho_integral(rho0,margi,margj,Z1,Z2,X1,X2,WIP,detJ) );
