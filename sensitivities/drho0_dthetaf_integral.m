function [drho0_dthetafi,drho0_dthetafj] = drho0_dthetaf_integral(rho,rho0,margi,margj,Z1,Z2,X1,X2,WIP,detJ)

Xi = X1; Zi = Z1;
Xj = X2; Zj = Z2;

F = (Xi-margi(2))/margi(3) .* (Xj-margj(2))/margj(3);

PHI2 = 1/(2*pi*sqrt(1-rho0^2)) * exp( -1/(2*(1-rho0^2)) * ( Zi.^2 - 2*rho0*Zi.*Zj + Zj.^2) );

dZi_dthetafi = dZ_dthetaf(Xi,margi);
dZj_dthetafj = dZ_dthetaf(Xj,margj);

dXi_dZi = dX_dZ(Xi,Zi,margi);
dXj_dZj = dX_dZ(Xj,Zj,margj);

dF_dXi = 1/margi(3)*(Xj-margj(2))/margj(3);
dF_dXj = (Xi-margi(2))/margi(3)/margj(3);

dF_dthetafi.mu    = -1/margi(3)*(Xj-margj(2))/margj(3);
dF_dthetafi.sigma = -(Xi-margi(2))/margi(3)^2.*(Xj-margj(2))/margj(3);
[dmui_dthetafi,dsigmai_dthetafi] = dmusigma_dp(margi);

dF_dthetafj.mu    = -(Xi-margi(2))/margi(3)/margj(3);
dF_dthetafj.sigma = -(Xi-margi(2))/margi(3).*(Xj-margj(2))/margj(3)^2;
[dmui_dthetafj,dsigmai_dthetafj] = dmusigma_dp(margj);

dF_dthetafi.mu    = -dF_dXi .* dXi_dZi .* dZi_dthetafi.mu + dF_dthetafi.mu;
dF_dthetafi.sigma = -dF_dXi .* dXi_dZi .* dZi_dthetafi.sigma + dF_dthetafi.sigma;
switch margi(1)
   case 8
      dF_dthetafi.p1 = dF_dthetafi.mu * dmui_dthetafi.p1 + dF_dthetafi.sigma * dsigmai_dthetafi.p1;
      dF_dthetafi.p2 = 0; dF_dthetafi.p3 = 0; dF_dthetafi.p4 = 0;
   case { 1, 2, 3, 4, 5, 6, 11, 12, 13, 15, 16 }
      dF_dthetafi.p1 = dF_dthetafi.mu * dmui_dthetafi.p1 + dF_dthetafi.sigma * dsigmai_dthetafi.p1;
      dF_dthetafi.p2 = dF_dthetafi.mu * dmui_dthetafi.p2 + dF_dthetafi.sigma * dsigmai_dthetafi.p2;
      dF_dthetafi.p3 = 0; dF_dthetafi.p4 = 0;
   case 14
      dF_dthetafi.p1 = dF_dthetafi.mu * dmui_dthetafi.p1 + dF_dthetafi.sigma * dsigmai_dthetafi.p1;
      dF_dthetafi.p2 = dF_dthetafi.mu * dmui_dthetafi.p2 + dF_dthetafi.sigma * dsigmai_dthetafi.p2;
      dF_dthetafi.p3 = dF_dthetafi.mu * dmui_dthetafi.p3 + dF_dthetafi.sigma * dsigmai_dthetafi.p3;
      dF_dthetafi.p4 = 0;
   case { 7, 51 }
      dF_dthetafi.p1 = dF_dthetafi.mu * dmui_dthetafi.p1 + dF_dthetafi.sigma * dsigmai_dthetafi.p1;
      dF_dthetafi.p2 = dF_dthetafi.mu * dmui_dthetafi.p2 + dF_dthetafi.sigma * dsigmai_dthetafi.p2;
      dF_dthetafi.p3 = dF_dthetafi.mu * dmui_dthetafi.p3 + dF_dthetafi.sigma * dsigmai_dthetafi.p3;
      dF_dthetafi.p4 = dF_dthetafi.mu * dmui_dthetafi.p4 + dF_dthetafi.sigma * dsigmai_dthetafi.p4;
end

dF_dthetafj.mu    = -dF_dXj .* dXj_dZj .* dZj_dthetafj.mu + dF_dthetafj.mu;
dF_dthetafj.sigma = -dF_dXj .* dXj_dZj .* dZj_dthetafj.sigma + dF_dthetafj.sigma;
switch margj(1)
   case 8
      dF_dthetafj.p1 = dF_dthetafj.mu * dmui_dthetafj.p1 + dF_dthetafj.sigma * dsigmai_dthetafj.p1;
      dF_dthetafj.p2 = 0; dF_dthetafj.p3 = 0; dF_dthetafj.p4 = 0;
   case { 1, 2, 3, 4, 5, 6, 11, 12, 13, 15, 16 }
      dF_dthetafj.p1 = dF_dthetafj.mu * dmui_dthetafj.p1 + dF_dthetafj.sigma * dsigmai_dthetafj.p1;
      dF_dthetafj.p2 = dF_dthetafj.mu * dmui_dthetafj.p2 + dF_dthetafj.sigma * dsigmai_dthetafj.p2;
      dF_dthetafj.p3 = 0; dF_dthetafj.p4 = 0;
   case 14
      dF_dthetafj.p1 = dF_dthetafj.mu * dmui_dthetafj.p1 + dF_dthetafj.sigma * dsigmai_dthetafj.p1;
      dF_dthetafj.p2 = dF_dthetafj.mu * dmui_dthetafj.p2 + dF_dthetafj.sigma * dsigmai_dthetafj.p2;
      dF_dthetafj.p3 = dF_dthetafj.mu * dmui_dthetafj.p3 + dF_dthetafj.sigma * dsigmai_dthetafj.p3;
      dF_dthetafj.p4 = 0;
   case { 7, 51 }
      dF_dthetafj.p1 = dF_dthetafj.mu * dmui_dthetafj.p1 + dF_dthetafj.sigma * dsigmai_dthetafj.p1;
      dF_dthetafj.p2 = dF_dthetafj.mu * dmui_dthetafj.p2 + dF_dthetafj.sigma * dsigmai_dthetafj.p2;
      dF_dthetafj.p3 = dF_dthetafj.mu * dmui_dthetafj.p3 + dF_dthetafj.sigma * dsigmai_dthetafj.p3;
      dF_dthetafj.p4 = dF_dthetafj.mu * dmui_dthetafj.p4 + dF_dthetafj.sigma * dsigmai_dthetafj.p4;
end

dPHI2_dZi   = -1./2./pi./(1-rho0.^2).^(1./2)./(2-2.*rho0.^2).*(2.*Zi-2.*rho0.*Zj).*exp(-1./(2-2.*rho0.^2).*(Zi.^2-2.*rho0.*Zi.*Zj+Zj.^2));
dPHI2_dZj   = -1./2./pi./(1-rho0.^2).^(1./2)./(2-2.*rho0.^2).*(-2.*rho0.*Zi+2.*Zj).*exp(-1./(2-2.*rho0.^2).*(Zi.^2-2.*rho0.*Zi.*Zj+Zj.^2));
dPHI2_drho0 = 1./2./pi./(1-rho0.^2).^(3./2).*exp(-1./(2-2.*rho0.^2).*(Zi.^2-2.*rho0.*Zi.*Zj+Zj.^2)).*rho0+1./2./pi./(1-rho0.^2).^(1./2).*(-4./(2-2.*rho0.^2).^2.*(Zi.^2-2.*rho0.*Zi.*Zj+Zj.^2).*rho0+2./(2-2.*rho0.^2).*Zi.*Zj).*exp(-1./(2-2.*rho0.^2).*(Zi.^2-2.*rho0.*Zi.*Zj+Zj.^2)); % Validé

drho0_dthetafi.mu    = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.mu,PHI2,F,dPHI2_dZi,dZi_dthetafi.mu,dPHI2_drho0,detJ,WIP);
drho0_dthetafi.sigma = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.sigma,PHI2,F,dPHI2_dZi,dZi_dthetafi.sigma,dPHI2_drho0,detJ,WIP);
switch margi(1)
   case 8
      drho0_dthetafi.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p1,PHI2,F,dPHI2_dZi,dZi_dthetafi.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p2 = 0; drho0_dthetafi.p3 = 0; drho0_dthetafi.p4 = 0;
   case { 1, 2, 3, 4, 5, 6, 11, 12, 13, 15, 16 }
      drho0_dthetafi.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p1,PHI2,F,dPHI2_dZi,dZi_dthetafi.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p2 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p2,PHI2,F,dPHI2_dZi,dZi_dthetafi.p2,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p3 = 0; drho0_dthetafi.p4 = 0;
   case 14
      drho0_dthetafi.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p1,PHI2,F,dPHI2_dZi,dZi_dthetafi.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p2 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p2,PHI2,F,dPHI2_dZi,dZi_dthetafi.p2,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p3 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p3,PHI2,F,dPHI2_dZi,dZi_dthetafi.p3,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p4 = 0;
   case { 7, 51 }
      drho0_dthetafi.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p1,PHI2,F,dPHI2_dZi,dZi_dthetafi.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p2 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p2,PHI2,F,dPHI2_dZi,dZi_dthetafi.p2,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p3 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p3,PHI2,F,dPHI2_dZi,dZi_dthetafi.p3,dPHI2_drho0,detJ,WIP);
      drho0_dthetafi.p4 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafi.p4,PHI2,F,dPHI2_dZi,dZi_dthetafi.p4,dPHI2_drho0,detJ,WIP);
end

drho0_dthetafj.mu    = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.mu,PHI2,F,dPHI2_dZj,dZj_dthetafj.mu,dPHI2_drho0,detJ,WIP);
drho0_dthetafj.sigma = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.sigma,PHI2,F,dPHI2_dZj,dZj_dthetafj.sigma,dPHI2_drho0,detJ,WIP);
switch margj(1)
   case 8
      drho0_dthetafj.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p1,PHI2,F,dPHI2_dZj,dZj_dthetafj.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p2 = 0; drho0_dthetafj.p3 = 0; drho0_dthetafj.p4 = 0;
   case { 1, 2, 3, 4, 5, 6, 11, 12, 13, 15, 16 }
      drho0_dthetafj.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p1,PHI2,F,dPHI2_dZj,dZj_dthetafj.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p2 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p2,PHI2,F,dPHI2_dZj,dZj_dthetafj.p2,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p3 = 0; drho0_dthetafj.p4 = 0;
   case 14
      drho0_dthetafj.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p1,PHI2,F,dPHI2_dZj,dZj_dthetafj.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p2 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p2,PHI2,F,dPHI2_dZj,dZj_dthetafj.p2,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p3 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p3,PHI2,F,dPHI2_dZj,dZj_dthetafj.p3,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p4 = 0;
   case { 7, 51 }
      drho0_dthetafj.p1 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p1,PHI2,F,dPHI2_dZj,dZj_dthetafj.p1,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p2 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p2,PHI2,F,dPHI2_dZj,dZj_dthetafj.p2,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p3 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p3,PHI2,F,dPHI2_dZj,dZj_dthetafj.p3,dPHI2_drho0,detJ,WIP);
      drho0_dthetafj.p4 = fzero('betadrho0_dthetaf',0,optimset('fzero'),dF_dthetafj.p4,PHI2,F,dPHI2_dZj,dZj_dthetafj.p4,dPHI2_drho0,detJ,WIP);
end