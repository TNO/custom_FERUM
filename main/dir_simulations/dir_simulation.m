function [ dirsimulationresults, probdata ] = dir_simulation(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

%     Finite Element Reliability Using Matlab, FERUM, Version 4.0, 2009 
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     A copy of the GNU General Public License is found in the file
%     <gpl-3.0.txt> following this collection of FERUM program files.
%     This license can be also found at: <http://www.gnu.org/licenses/>.    
%
%     For more information on FERUM, visit: <http://www.ifma.fr/FERUM/>


global nfun
nfun = 0;

% Extract model data
marg = probdata.marg;
R = probdata.correlation;
transf_type = probdata.transf_type;

if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end

block_size     = analysisopt.block_size;

rand_generator = analysisopt.rand_generator;
stdv1          = analysisopt.stdv_sim;
num_sim        = analysisopt.num_sim;
target_cov     = analysisopt.target_cov;

dir_flag       = analysisopt.dir_flag;


% Find number of random variables
nrv = size(marg,1);
% Set point for directional simulations
point = zeros(nrv,1); 

rho = analysisopt.rho;
if isfield(analysisopt,'tolx')
   tolx = analysisopt.tolx;
else
   tolx = 1e-5;
end
if isfield(analysisopt,'keep_a') & analysisopt.keep_a == 1
   keep_a = 1;
   % Initialization for storage of a values which gives axes along which simulations are carried out
   Alla = zeros(nrv,num_sim);
else
   keep_a = 0;
end
if isfield(analysisopt,'keep_r') & analysisopt.keep_r == 1
   keep_r = 1;
   % Initialization for storage of r values for which g(r) = 0
   Allr = zeros(1,num_sim);
else
   keep_r = 0;
end


% Modify correlation matrix and perform Cholesky decomposition
if ~isfield(probdata,'Lo')
    
   if transf_type == 3
    
      % Compute corrected correlation coefficients
      switch probdata.Ro_method
         case 0
            Ro = mod_corr( marg, R );
         case 1
            [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , 0);
      end     
      probdata.Ro = Ro;
   
      % Cholesky decomposition
      [Lo,ierr] = my_chol(Ro); if  ierr>0, return, end
      probdata.Lo = Lo;
   
      iLo = inv(Lo);
      probdata.iLo = iLo;

   end

end


% Check for large failure probabilities
if sigmafun(0,ones(nrv,1),rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield) <= 0
   
   dirsimulationresults.pf           = 1;
   dirsimulationresults.cov_pf       = nan;
   dirsimulationresults.beta         = -inf;
   dirsimulationresults.npts         = 1;
   if keep_a
      dirsimulationresults.Alla      = ones(nrv,1);
   end
   if keep_r
      dirsimulationresults.Allr      = 0;
   end
   dirsimulationresults.nfun         = nfun;
    
   return
    
end


% Establish covariance matrix, its Cholesky decomposition, and its inverse
covariance = stdv1^2 * eye(nrv);
chol_covariance = chol(covariance);
inv_covariance = inv(covariance);

% Initializations
sum_q = 0;
sum_q_squared = 0;
cov_of_q_bar = nan*ones(1,num_sim);

% Pre-compute some factors to minimize computations inside simulation loop
factor1 = 1 / ( (2*pi)^(nrv/2) );
factor2 = 1 / ( (2*pi)^(nrv/2) * sqrt(det(covariance)) );

k = 0;
cov_of_q_bar(1) = 1.0;
percent_done = 0;

switch dir_flag
   % Partition the unit hypersphere into regions of equal area - Directions for directional simulations
   % given by the centers of these regions
   case 'det'
      Allu = point*ones(1,num_sim) + chol_covariance * eq_point_set(nrv-1,num_sim);
end


while k < num_sim
   
   block_size = min( block_size, num_sim-k );
   k = k + block_size;
   
   % Do the simulation (create array of random numbers)
   switch dir_flag

      case 'random'

         switch rand_generator
            case 0
               allu = point*ones(1,block_size) + chol_covariance * randn(nrv,block_size);
            otherwise
%                allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(twister(nrv,block_size));
               allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(rand(nrv,block_size));
         end

      case 'det'

         allu = Allu(:,(k-block_size+1):k);

   end
   
   % Find normalized u - Unit vectors of the directional axes
   alla = allu ./ ( ones(nrv,1) * dot(allu,allu).^0.5 );
   if keep_a
      Alla(:,(k-block_size+1):k) = alla;
   end

   % fzero search for g = 0 in [0 rho] interval
   [ allr, allfval ] = my_fzero_sigmafun_vectorized([0 rho+tolx*10]'*ones(1,block_size),tolx,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
   if keep_r
      Allr((k-block_size+1):k) = allr;
   end
  
   % Compute values of joint distributions at the u-point
   allphi = factor1 * exp( -0.5 * dot(allu,allu) );
   allh   = factor2 * exp( -0.5 * dot ( (allu-point*ones(1,block_size)) , inv_covariance * (allu-point*ones(1,block_size)) ) );

   % Update sums
   allq = ( 1 - ferum_cdf(8,allr.^2,[nrv, NaN, NaN, 1]) ) .* allphi ./ allh;
   sum_q = sum_q + sum(allq);
   sum_q_squared = sum_q_squared + sum(allq.^2);
   
   % Compute coefficient of variation (of pf) for "block_size" simulations
   if sum_q > 0
      q_bar(k) = 1/k * sum_q;
      variance_of_q_bar = 1/k * ( 1/k * sum_q_squared - (1/k*sum_q)^2 );
      cov_of_q_bar(k) = sqrt(variance_of_q_bar) / q_bar(k);
      if cov_of_q_bar(k) == 0
         cov_of_q_bar(k) = 1.0;
      end
   else
      q_bar(k) = 0;
      cov_of_q_bar(k) = 1.0;
   end

   if floor( k/num_sim * 20 ) > percent_done
      if echo_flag
         percent_done = floor( k/num_sim * 20 );
         fprintf(1,'%d%% complete\n',percent_done*5)
      end
   end
      
end


if sum_q > 0
   
   % Compute probability of failure and reliability index
   pf   = q_bar(k);
   cov  = cov_of_q_bar(k);
   beta = -inv_norm_cdf(pf);
   
   if echo_flag
       
      figure, hold on
      if ( nrv <= 3 ) & ( num_sim <= 10000 ) & keep_a & keep_r
         X = (ones(nrv,1)*Allr) .* Alla;
         Inan = find( abs ( Allr ) > ( rho - 10*tolx ) );
         X(:,Inan) = nan*ones(nrv,length(Inan));
         [ rmin, Imin ] = min( Allr ); % X(:,Imin)
         if nrv == 3
            plot3(X(1,:),X(2,:),X(3,:),'b.')
            hold on
            plot3(X(1,Imin),X(2,Imin),X(3,Imin),'go','MarkerSize',8,'MarkerFaceColor',[0 1 0])
            plot3([0 X(1,Imin)],[0 X(2,Imin)],[0 X(3,Imin)],'g-','MarkerSize',8,'MarkerFaceColor',[0 1 0])
            plot3(0,0,0,'ks','MarkerSize',8,'MarkerFaceColor',[0 0 0])
            xlabel(probdata.name(1),'FontSize',18); ylabel(probdata.name(2),'FontSize',18); zlabel(probdata.name(3),'FontSize',18);
            set(gca,'FontSize',18);
            axis equal
            grid on
            box on
            title('U-space','FontSize',18)
            rotate3d on
         elseif nrv == 2
            plot(X(1,:),X(2,:),'b.')
            hold on
            plot(X(1,Imin),X(2,Imin),'go','MarkerSize',8,'MarkerFaceColor',[0 1 0])
            plot([0 X(1,Imin)],[0 X(2,Imin)],'g-','MarkerSize',8,'MarkerFaceColor',[0 1 0])
            plot(0,0,'ks','MarkerSize',8,'MarkerFaceColor',[0 0 0])
            xlabel(probdata.name(1),'FontSize',18); ylabel(probdata.name(2),'FontSize',18);
            set(gca,'FontSize',18);
            axis equal
            grid on
            box on
            title('U-space','FontSize',18)
         end
      end
      
   end
   
else
   
   pf   = 0;
   cov  = 0;
   beta = 0;
   
end

dirsimulationresults.pf      = pf;
dirsimulationresults.beta    = beta;
dirsimulationresults.ndir    = k;
switch dir_flag
   case 'random'
      dirsimulationresults.cov_pf  = cov;
end
if keep_a
   dirsimulationresults.Alla = Alla;
end
if keep_r
   dirsimulationresults.Allr = Allr;
end
dirsimulationresults.nfun    = nfun;