function [ distributionresults, probdata ] = distribution_analysis(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

% Find number of random variables
nrv = size(marg,1);

% Set point for crude Monte Carlo / importance sampling 
switch lower(analysisopt.sim_point)
   case 'origin'
      point = zeros(nrv,1);
   case 'dspt'
      point = analysisopt.formresults.dsptu;
end

ALLX = zeros(nrv,num_sim);
if isfield(gfundata(lsf),'ng')
   ng = gfundata(lsf).ng;
else
   ng = 1;
end
ALLG = zeros(ng,num_sim);

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
chol_covariance = chol(covariance);

k = 0;
percent_done = 0;

while k < num_sim
   
   block_size = min( block_size, num_sim-k );
   k = k + block_size;
   
   % Do the simulation (create array of random numbers)
   switch rand_generator
      case 0
         allu = point*ones(1,block_size) + chol_covariance * randn(nrv,block_size);
       otherwise
           %replacement
           allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(rand(nrv,block_size));
           
           %original
           %allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(twister(nrv,block_size));
   end
   
   % Transform into original space
   allx = u_to_x(allu,probdata);
   ALLX(:,(k-block_size+1):k) = allx;

   % Evaluate limit-state function
   [ allg, dummy ] = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
   ALLG(:,(k-block_size+1):k) = allg;
   
   if echo_flag
      if floor( k/num_sim * 20 ) > percent_done
         percent_done = floor( k/num_sim * 20 );
         fprintf(1,'%d%% complete\n',percent_done*5)
      end
   end
      
end

% c.f. whisto.m from wafo toolbox - Call whisto(xi,N,0,1)
%
% N = approximate number of bins wanted ( default depending on size(ALLX,2) ) 
N = ceil(4*sqrt(sqrt(size(ALLX,2))));
% odd = placement of bins (0 or 1) (default 0)
odd = 0;
% scale = argument for scaling (default 0)
%         scale = 1 yields the area 1 under the histogram
scale = 1;

if echo_flag
   
   while 1

      close all
      
      for i = 1:nrv

         xi = ALLX(i,:);
         [mn, mx, d] = plot_histo(i,xi,N,odd,scale);

         xxref = linspace(mn,mx,200);
         pdfref = ferum_pdf(marg(i,1),xxref, marg(i,5:8));
         plot(xxref,pdfref,'r-')
         set(gca,'FontSize',18);

         xlabel(probdata.name(i),'FontSize',18);
         ylabel('Probability Density','FontSize',18);

      end

      for ig=1:size(ALLG,1)

         xr = ALLG(ig,:);
         [mn, mx, d] = plot_histo(nrv+ig,xr,N,odd,scale);

         xlabel(sprintf('Output g(%d)',ig));
         ylabel('Probability Density');

      end

      fprintf(1,'\nPick a new value for approximated number of bins wanted\n( press enter to keep current value %d )', N)
      N = input('\nN = ');
      if isempty(N), break, end

   end
   
end


if echo_flag
   disp([' '])
   disp(['..............................................................................................'])
   disp([' '])
end
for i = 1:nrv
   xi = ALLX(i,:);
   Xmean(i,1) = mean(xi); Xstdv(i,1) = std(xi);
   if echo_flag
      fprintf(1, '%6s :  mu = %12.4e (ref = %12.4e)  -  stdv = %12.4e (ref = %12.4e)\n', ...
         char(probdata.name(i)), mean(xi), marg(i,2), std(xi), marg(i,3));
   end
end
fprintf('\n');
for ig=1:size(ALLG,1)
   xr = ALLG(ig,:);
   Gmean = mean(xr); Gstdv = std(xr);
   if echo_flag
      fprintf(1, 'Output g(%d) :  mu = %12.4e                       -  stdv = %12.4e\n', ...
              ig, mean(xr), std(xr));
   end
end


distributionresults.Xmean = Xmean;
distributionresults.Xstdv = Xstdv;
distributionresults.Gmean = Gmean;
distributionresults.Gstdv = Gstdv;
distributionresults.ALLX  = ALLX;
distributionresults.ALLG  = ALLG;
distributionresults.nfun  = k;
