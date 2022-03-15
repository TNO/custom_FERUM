function [ dbeta_drho, dpf_drho, dbeta_dthetaf, dpf_dthetaf, dbeta_dthetaf_approx, dbeta_dthetaf_ratio ] = sensitivities_wrt_thetaf( beta, alpha, x, u, marg, Ro, Lo, iLo, dRo_drho, dRo_dthetafi, dRo_dthetafj )

nrv = size(marg,1);

z = Lo * u;


% Beta and pf sensitivities w.r.t. correlation coefficients

dbeta_drho = zeros(nrv);
dpf_drho   = zeros(nrv);

if ~isempty(dRo_drho)

   for i = 1 : nrv
      for j = 1 : (i-1)
         dRo_drhoij = zeros(nrv);
         dRo_drhoij(i,j) = dRo_drho(i,j); dRo_drhoij(j,i) = dRo_drho(j,i);
         [Lo,dLo_drhoij,ierr] = my_dchol(Ro,dRo_drhoij);
         diLo_drhoij = -iLo*dLo_drhoij*iLo;
         dbeta_drho(i,j) = alpha' * ( diLo_drhoij * z );
      end
   end
   dbeta_drho = dbeta_drho + tril(dbeta_drho,-1)';

   % pf sensitivities w.r.t. correlation coefficients
   dpf_drho = -normpdf(beta) * dbeta_drho;

end

% Beta sensitivities w.r.t. mean, standard deviation and distribution parameters

dbeta_dthetaf        = zeros(nrv,6);
dbeta_dthetaf_approx = zeros(nrv,6);
dbeta_dthetaf_ratio  = zeros(nrv,6);
dpf_dthetaf          = zeros(nrv,6);

% Compute the partial derivatives of the transformation w.r.t. distribution parameters
J_z_thetaf = jacobian_z_thetaf(x,marg);

for i = 1 : nrv
    
   % Beta sensitivities w.r.t. mean
   dbeta_dthetaf_approx(i,1) = alpha' * iLo * J_z_thetaf(:,6*(i-1)+1);
   if ~isempty(dRo_dthetafi)
      dRo_dthetaf = zeros(nrv);
      dRo_dthetaf(i,:) = dRo_dthetafi.mu(i,:);
      dRo_dthetaf(:,i) = dRo_dthetafj.mu(:,i);
      dRo_dthetaf = dRo_dthetaf + tril(dRo_dthetaf,-1)';
   
      [Lo,dLo_dthetaf,ierr] = my_dchol(Ro,dRo_dthetaf);
      diLo_dthetaf = -iLo*dLo_dthetaf*iLo;
   
      dbeta_dthetaf(i,1) = alpha' * ( diLo_dthetaf * z + iLo * J_z_thetaf(:,6*(i-1)+1) );
      if abs(dbeta_dthetaf(i,1)) < eps
         dbeta_dthetaf_ratio(i,1) = nan;
      else
         dbeta_dthetaf_ratio(i,1) = ( alpha' * diLo_dthetaf * z ) / ( alpha' * iLo * J_z_thetaf(:,6*(i-1)+1) );
      end
   end
   
   % Beta sensitivities w.r.t. standard deviation
   dbeta_dthetaf_approx(i,2) = alpha' * iLo * J_z_thetaf(:,6*(i-1)+2);
   if ~isempty(dRo_dthetafi)
      dRo_dthetaf = zeros(nrv);
      dRo_dthetaf(i,:) = dRo_dthetafi.sigma(i,:);
      dRo_dthetaf(:,i) = dRo_dthetafj.sigma(:,i);
      dRo_dthetaf = dRo_dthetaf + tril(dRo_dthetaf,-1)';

      [Lo,dLo_dthetaf,ierr] = my_dchol(Ro,dRo_dthetaf);
      diLo_dthetaf = -iLo*dLo_dthetaf*iLo;

      dbeta_dthetaf(i,2) = alpha' * ( diLo_dthetaf * z + iLo * J_z_thetaf(:,6*(i-1)+2) );
      if abs(dbeta_dthetaf(i,2)) < eps
         dbeta_dthetaf_ratio(i,2) = nan;
      else
         dbeta_dthetaf_ratio(i,2) = ( alpha' * diLo_dthetaf * z ) / ( alpha' * iLo * J_z_thetaf(:,6*(i-1)+2) );
      end
   end
   
   % Beta sensitivities w.r.t. p1 parameter
   dbeta_dthetaf_approx(i,3) = alpha' * iLo * J_z_thetaf(:,6*(i-1)+3);
   if ~isempty(dRo_dthetafi)
      dRo_dthetaf = zeros(nrv);
      dRo_dthetaf(i,:) = dRo_dthetafi.p1(i,:);
      dRo_dthetaf(:,i) = dRo_dthetafj.p1(:,i);
      dRo_dthetaf = dRo_dthetaf + tril(dRo_dthetaf,-1)';

      [Lo,dLo_dthetaf,ierr] = my_dchol(Ro,dRo_dthetaf);
      diLo_dthetaf = -iLo*dLo_dthetaf*iLo;

      dbeta_dthetaf(i,3) = alpha' * ( diLo_dthetaf * z + iLo * J_z_thetaf(:,6*(i-1)+3) );
      if abs(dbeta_dthetaf(i,3)) < eps
         dbeta_dthetaf_ratio(i,3) = nan;
      else
         dbeta_dthetaf_ratio(i,3) = ( alpha' * diLo_dthetaf * z ) / ( alpha' * iLo * J_z_thetaf(:,6*(i-1)+3) );
      end
   end

   if marg(i,1) ~= 8 % Beta sensitivities w.r.t. p2 parameter
   
   dbeta_dthetaf_approx(i,4) = alpha' * iLo * J_z_thetaf(:,6*(i-1)+4);
   if ~isempty(dRo_dthetafi)
      dRo_dthetaf = zeros(nrv);
      dRo_dthetaf(i,:) = dRo_dthetafi.p2(i,:);
      dRo_dthetaf(:,i) = dRo_dthetafj.p2(:,i);
      dRo_dthetaf = dRo_dthetaf + tril(dRo_dthetaf,-1)';

      [Lo,dLo_dthetaf,ierr] = my_dchol(Ro,dRo_dthetaf);
      diLo_dthetaf = -iLo*dLo_dthetaf*iLo;

      dbeta_dthetaf(i,4) = alpha' * ( diLo_dthetaf * z + iLo * J_z_thetaf(:,6*(i-1)+4) );
      if abs(dbeta_dthetaf(i,4)) < eps
         dbeta_dthetaf_ratio(i,4) = nan;
      else
         dbeta_dthetaf_ratio(i,4) = ( alpha' * diLo_dthetaf * z ) / ( alpha' * iLo * J_z_thetaf(:,6*(i-1)+4) );
      end
   end
   
   end % Beta sensitivities w.r.t. p2 parameter

   if marg(i,1) == 14 | marg(i,1) == 7 % Beta sensitivities w.r.t. p3 parameter
       
   dbeta_dthetaf_approx(i,5) = alpha' * iLo * J_z_thetaf(:,6*(i-1)+5);
   if ~isempty(dRo_dthetafi)
      dRo_dthetaf = zeros(nrv);
      dRo_dthetaf(i,:) = dRo_dthetafi.p3(i,:);
      dRo_dthetaf(:,i) = dRo_dthetafj.p3(:,i);
      dRo_dthetaf = dRo_dthetaf + tril(dRo_dthetaf,-1)';

      [Lo,dLo_dthetaf,ierr] = my_dchol(Ro,dRo_dthetaf);
      diLo_dthetaf = -iLo*dLo_dthetaf*iLo;

      dbeta_dthetaf(i,5) = alpha' * ( diLo_dthetaf * z + iLo * J_z_thetaf(:,6*(i-1)+5) );
      if abs(dbeta_dthetaf(i,5)) < eps
         dbeta_dthetaf_ratio(i,5) = nan;
      else
         dbeta_dthetaf_ratio(i,5) = ( alpha' * diLo_dthetaf * z ) / ( alpha' * iLo * J_z_thetaf(:,6*(i-1)+5) );
      end
   end

   end % Beta sensitivities w.r.t. p3 parameter

   if marg(i,1) == 7 % Beta sensitivities w.r.t. p4 parameter
       
   dbeta_dthetaf_approx(i,6) = alpha' * iLo * J_z_thetaf(:,6*(i-1)+6);
   if ~isempty(dRo_dthetafi)
      dRo_dthetaf = zeros(nrv);
      dRo_dthetaf(i,:) = dRo_dthetafi.p4(i,:);
      dRo_dthetaf(:,i) = dRo_dthetafj.p4(:,i);
      dRo_dthetaf = dRo_dthetaf + tril(dRo_dthetaf,-1)';

      [Lo,dLo_dthetaf,ierr] = my_dchol(Ro,dRo_dthetaf);
      diLo_dthetaf = -iLo*dLo_dthetaf*iLo;

      dbeta_dthetaf(i,6) = alpha' * ( diLo_dthetaf * z + iLo * J_z_thetaf(:,6*(i-1)+6) );
      if abs(dbeta_dthetaf(i,6)) < eps
         dbeta_dthetaf_ratio(i,6) = nan;
      else
         dbeta_dthetaf_ratio(i,6) = ( alpha' * diLo_dthetaf * z ) / ( alpha' * iLo * J_z_thetaf(:,6*(i-1)+6) );
      end
   end
   
   end % Beta sensitivities w.r.t. p4 parameter
                                           
end

if isempty(dRo_dthetafi)
   % Approximated sensitivities
   dbeta_dthetaf = dbeta_dthetaf_approx;
end

% pf sensitivities w.r.t. mean, standard deviation and distribution parameters
dpf_dthetaf = -normpdf(beta) * dbeta_dthetaf;
