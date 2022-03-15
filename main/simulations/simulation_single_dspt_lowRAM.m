function [ simulationresults, probdata ] = simulation_single_dspt_lowRAM(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

block_size = analysisopt.block_size;

rand_generator = analysisopt.rand_generator;
stdv1 = analysisopt.stdv_sim;
num_sim = analysisopt.num_sim;
target_cov = analysisopt.target_cov;

% Find number of random variables
nrv = size(marg,1);

% Set point for crude Monte Carlo / importance sampling 
switch lower(analysisopt.sim_point)
   case 'origin'
      point = zeros(nrv,1);
   case 'dspt'
      point = analysisopt.formresults.dsptu;
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

% Establish covariance matrix, its Cholesky decomposition, and its inverse
covariance = stdv1^2 * eye(nrv);
chol_covariance = stdv1 * eye(nrv);  % chol_covariance = chol(covariance);
inv_covariance = 1/stdv1^2 * eye(nrv); % inv_covariance = inv(covariance);

% Initializations
sum_q = 0;
sum_q_squared = 0;

% Pre-compute some factors to minimize computations inside simulation loop
factor1_over_factor2 = stdv1^nrv;

% Initializations
k = 0;
cov_of_q_bar(1) = 1.0;
percent_done = 0;

while k < num_sim
   
   block_size = min( block_size, num_sim-k );
   k = k + block_size;
   
   % Do the simulation (create array of random numbers)
   switch rand_generator
      case 0
         allu = point*ones(1,block_size) + chol_covariance * randn(nrv,block_size);
      otherwise
         allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(twister(nrv,block_size));
   end
   
   % Transform into original space
   allx = u_to_x(allu,probdata);
   
   % Evaluate limit-state function
   [ allg, dummy ] = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
  
   % Collect result of sampling: if g < 0 , I = 1 , else I = 0
   allI = zeros(1,block_size);
   Ifail = find( allg < 0 );
   allI(Ifail) = 1;
   
   % Update sums
   allq = allI * factor1_over_factor2 .* exp( -0.5 * dot(allu,allu) + ...
                                              0.5 * dot ( (allu-point*ones(1,block_size)) , inv_covariance * (allu-point*ones(1,block_size)) ) );
   sum_q = sum_q + sum(allq);
   sum_q_squared = sum_q_squared + sum(allq.^2);
   
   % Compute coefficient of variation (of pf) for block_size samples
   if sum_q > 0
      q_bar = 1/k * sum_q;
      variance_of_q_bar = 1/k * ( 1/k * sum_q_squared - (1/k*sum_q)^2 );
      cov_of_q_bar = sqrt(variance_of_q_bar) / q_bar;
      if cov_of_q_bar == 0
         cov_of_q_bar = 1.0;
      end
   else
      q_bar = 0;
      cov_of_q_bar = 1.0;
   end

   if floor( k/num_sim * 20 ) > percent_done
      percent_done = floor( k/num_sim * 20 );
      if echo_flag
         fprintf(1,'%d%% complete\n',percent_done*5)
      end
   end
      
   if cov_of_q_bar <= target_cov, break, end
   
end


if sum_q > 0
   
   % Compute probability of failure and reliability index
   pf = q_bar;
   cov = cov_of_q_bar;
   beta = -inv_norm_cdf(pf);
      
else
   
   pf = 0;
   cov = 0;
   beta = 0;
   
end

simulationresults.pf           = pf;
simulationresults.cov_pf       = cov;
simulationresults.beta         = beta;
simulationresults.nfun         = k;