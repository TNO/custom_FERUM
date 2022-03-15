function [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , flag_sens)

% MOD_CORR_SOLVE    Modifies correlation matrix
%
%   [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , flag_sens )
%
%   Function 'mod_corr_solve' modifies the correlation matrix according to the assumption
%   of Nataf distributed random variables (bifolded integral equation solved numerically).


nrv = size(marg,1);


Ro = eye(size(R));


if flag_sens == 1
    
   dRo_drho = eye(size(R));
   
   dRo_dthetafi.mu    = zeros(size(R));
   dRo_dthetafi.sigma = zeros(size(R));
   dRo_dthetafi.p1    = zeros(size(R));
   dRo_dthetafi.p2    = zeros(size(R));
   dRo_dthetafi.p3    = zeros(size(R));
   dRo_dthetafi.p4    = zeros(size(R));

   dRo_dthetafj.mu    = zeros(size(R));
   dRo_dthetafj.sigma = zeros(size(R));
   dRo_dthetafj.p1    = zeros(size(R));
   dRo_dthetafj.p2    = zeros(size(R));
   dRo_dthetafj.p3    = zeros(size(R));
   dRo_dthetafj.p4    = zeros(size(R));

else
    
   dRo_drho = []; dRo_dthetafi = []; dRo_dthetafj = [];
    
end


for i = 1 : nrv       % Loop through all off-diagonal elements of the correlation matrix

   for j = 1:(i-1) 
                         
      rho = R(i,j);

      if rho ~= 0 | flag_sens
          
         margi = marg(i,:);
         margj = marg(j,:);
         
         zmax = 6;
         
         if abs(rho) > 0.9995
            nIP = 1024;
         elseif abs(rho) > 0.998
            nIP = 512;
         elseif abs(rho) > 0.992
            nIP = 256;
         elseif abs(rho) > 0.97
            nIP = 128;
         elseif abs(rho) > 0.9
            nIP = 64;
         else
            nIP = 32;
         end
        
         [Z1,Z2,X1,X2,WIP,detJ] = zi_and_xi_points(margi,margj,zmax,nIP);
         
      end
     
      if rho ~= 0
         rho0 = fminsearch('abs_of_betacorr',rho,optimset('fzero'),rho,margi,margj,Z1,Z2,X1,X2,WIP,detJ);
         if isnan(rho0) % In general, never happens!
             rrho = -0.99:0.01:0.99; ff = zeros(size(rrho));
             for mm = 1:length(rrho)
                ff(mm) = betacorr(rrho(mm),rho,margi,margj,Z1,Z2,X1,X2,WIP,detJ);
             end
             figure(1)
             plot(rrho,ff)
             xlabel('rho0');
             ylabel(sprintf('%.4f - rho(rho0)',rho));
             grid on
         end
      else
         rho0 = 0;
      end
     
      Ro(i,j) = rho0;
      
      if flag_sens == 1
             
         % Sensitivities of R0 w.r.t. correlation coefficients 
         drho0_drho = drho0_drho_integral(rho0,margi,margj,Z1,Z2,X1,X2,WIP,detJ);
         dRo_drho(i,j) = drho0_drho;
         
         % Sensitivities of R0 w.r.t. distribution parameters (Gauss integration scheme)
         [drho0_dthetafi_Gauss,drho0_dthetafj_Gauss] = drho0_dthetaf_integral(rho,rho0,margi,margj,Z1,Z2,X1,X2,WIP,detJ);
         dRo_dthetafi.mu(i,j)    = drho0_dthetafi_Gauss.mu;
         dRo_dthetafi.sigma(i,j) = drho0_dthetafi_Gauss.sigma;
         dRo_dthetafi.p1(i,j)    = drho0_dthetafi_Gauss.p1;
         dRo_dthetafi.p2(i,j)    = drho0_dthetafi_Gauss.p2;
         dRo_dthetafi.p3(i,j)    = drho0_dthetafi_Gauss.p3;
         dRo_dthetafi.p4(i,j)    = drho0_dthetafi_Gauss.p4;
         dRo_dthetafj.mu(i,j)    = drho0_dthetafj_Gauss.mu;
         dRo_dthetafj.sigma(i,j) = drho0_dthetafj_Gauss.sigma;
         dRo_dthetafj.p1(i,j)    = drho0_dthetafj_Gauss.p1;
         dRo_dthetafj.p2(i,j)    = drho0_dthetafj_Gauss.p2;
         dRo_dthetafj.p3(i,j)    = drho0_dthetafj_Gauss.p3;
         dRo_dthetafj.p4(i,j)    = drho0_dthetafj_Gauss.p4;
         
      end % if flag_sens
         
   end % for j = 1:(i-1)
   
end % for i = 1 : nrv


Ro = Ro + tril(Ro,-1)';


if flag_sens == 1
   dRo_drho = dRo_drho + tril(dRo_drho,-1)';
end