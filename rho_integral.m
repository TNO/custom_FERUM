function rho = rho_integral(rho0,margi,margj,Z1,Z2,X1,X2,WIP,detJ)

% computes the bi-folded rho-integral by 2D numerical integration

% Previous non-vectorized version
%rho = 0;
%for i = 1:nIP
%   for j = 1:nIP
%      phi2(i,j) = 1/(2*pi*sqrt(1-rho0^2)) * exp( -1/(2*(1-rho0^2)) * ( z1(i)^2 - 2*rho0*z1(i)*z2(j) + z2(j)^2) );
%      rho = rho + (x1(i)-margi(2))/margi(3) * (x2(j)-margj(2))/margj(3) * phi2(i,j) * detJ * wIP(i) * wIP(j);
%   end
%end

phi2 = 1/(2*pi*sqrt(1-rho0^2)) * exp( -1/(2*(1-rho0^2)) * ( Z1.^2 - 2*rho0*Z1.*Z2 + Z2.^2) );
rho = sum ( sum ( (X1-margi(2))/margi(3) .* (X2-margj(2))/margj(3) .* phi2 * detJ .* WIP ) );
