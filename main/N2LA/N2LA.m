function [ N2LAresults, probdata, rbdo_parameters] = N2LA(lsf,rbdo_parameters,rbdo_fundata,probdata,analysisopt,gfundata,femodel,randomfield)

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


stepbystep_figures = 0;

femodel.save_formresults = 0;   % Save formresults per RBDO_iter = 1 / Don't save formresults per RBDO_iter = 0
femodel.RBDO_iter        = 0;   % For the very first FORM evaluation...

global rbdo_nform
rbdo_nform = 0;

global rbdo_nfun
rbdo_nfun = 0;

p = length(rbdo_fundata.cost);
q = length(rbdo_fundata.constraint) + 1;

rbdo_parameters.problemsize = [p q];
rbdo_parameters.warmstart   = [];

p = 1; % p in the rbdo file is equal to one like 'one objective function' returned by fun: c_0 + c_f * p_f

max_iter = rbdo_parameters.max_iter;

gfundata(lsf).thetag0 = gfundata(lsf).thetag;
thetag = ones(size(gfundata(lsf).thetag));

gfundata(lsf).thetag = thetag;
rbdo_funvalues_x.initial.f0 = ones(q,1);
rbdo_funvalues_x.initial.c0 = ones(p,1);

nthetag = length(thetag);

[ c, f ] = fun(lsf,thetag,'func',rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield);

very_initial_cost = c;

rbdo_funvalues_x.initial.f0 = abs(f);
rbdo_funvalues_x.initial.c0 = c;


initial_cost = c./rbdo_funvalues_x.initial.c0;
phi = max(c./rbdo_funvalues_x.initial.c0);
psi = max(f./rbdo_funvalues_x.initial.f0);
rbdo_funvalues_x.c   = c./rbdo_funvalues_x.initial.c0;
rbdo_funvalues_x.phi = phi;
rbdo_funvalues_x.f   = f./rbdo_funvalues_x.initial.f0;
rbdo_funvalues_x.psi = psi;

H = [ thetag' thetag'.*gfundata(lsf).thetag0' ...
      rbdo_funvalues_x.c'      rbdo_funvalues_x.c'.*rbdo_funvalues_x.initial.c0' ...
      rbdo_funvalues_x.phi      max(rbdo_funvalues_x.c'.*rbdo_funvalues_x.initial.c0') ...
      rbdo_funvalues_x.f'      rbdo_funvalues_x.f'.*rbdo_funvalues_x.initial.f0' ...
      rbdo_funvalues_x.psi     max(rbdo_funvalues_x.f'.*rbdo_funvalues_x.initial.f0') ...
      nan ]';

i = 1;
termflag = 1;

while termflag > 0

   precision_test = 1;

   while precision_test  > 0

      femodel.RBDO_iter = i;

      [Hinner, stepsize_k, gfundata, rbdo_funvalues_x] = ph_quadprog(thetag,1,rbdo_parameters,rbdo_funvalues_x,rbdo_fundata,lsf,probdata,analysisopt,gfundata,femodel,randomfield);

      z         = Hinner(1:nthetag, 2);
      cz        = Hinner(2*nthetag+1:2*nthetag+p, 2);
      phiz      = Hinner(2*nthetag+2*p+1, 2);
      fz        = Hinner(2*nthetag+2*p+2*1+1:2*nthetag+2*p+2*1+q, 2);
      psiz      = Hinner(2*nthetag+2*p+2*1+2*q+1, 2);
      psiz_plus = max(0,psiz);

      rbdo_parameters.method;
      thetag = z;

      rbdo_funvalues_x.c        = cz;
      rbdo_funvalues_x.phi      = phiz;
      rbdo_funvalues_x.f        = fz;
      rbdo_funvalues_x.psi      = psiz;

      rbdo_parameters.warmstart = stepsize_k;

      H = [H [Hinner(:,2)]]; % save H.mat H;
      precision_test = -1;

   end

   i = i + 1;

   disp(['..............................................................................................'])
   disp('RBDO analysis is running, please wait... (Ctrl+C breaks)')
   disp([' '])
   fprintf('Iteration # %d /%d\n\n',i-1,max_iter);
   for a=1:nthetag
      fprintf('%s = %d\n',gfundata(lsf).thetagname{a},thetag(a).*gfundata(lsf).thetag0(a));
   end
   fprintf('\nWith normalized cost   = %d\n', cz.*rbdo_funvalues_x.initial.c0./very_initial_cost);
   fprintf('\nWith real cost         = %d\n', cz.*rbdo_funvalues_x.initial.c0);
   fprintf('\nWith reliability index = %d\n', - rbdo_funvalues_x.f(q,1).*rbdo_funvalues_x.initial.f0(q,1) + rbdo_parameters.target_beta);
   fprintf('\nMethod: %s\n',rbdo_parameters.method);
   disp(['..............................................................................................'])

   if stepbystep_figures

      close all
      figure

      subplot(2,4,1)
      if i == 2, clf, end
      hold on
      plot((0:size(H,2)-1)', H(1:nthetag,:)','.-');
      title('Normalized design variables (\theta_g) / Iterations')
      xlabel('Iterations');
      ylabel('Design variables (\theta_g)');
      legend(char(gfundata(1).thetagname))

      subplot(2,4,5)
      if i == 2, clf, end
      hold on
      plot((0:size(H,2)-1)', H(nthetag+1:2*nthetag,:)','.-');
      title('Design variables (\theta_g) / Iterations')
      xlabel('Iterations');
      ylabel('Design variables (\theta_g)');
      legend(char(gfundata(1).thetagname))

      subplot(2,4,2)
      if i == 2, clf, end
      plot((0:size(H,2)-1)', H(2*nthetag+p+1:2*nthetag+2*p,:)'./very_initial_cost, 'r.-' );
      title('Normalized cost / Iterations')
      xlabel('Iterations');
      ylabel('Normalized cost');

      subplot(2,4,6)
      if i == 2, clf, end
      plot((0:size(H,2)-1)', H(2*nthetag+p+1:2*nthetag+2*p,:)', 'r.-' );
      title('Cost / Iterations')
      xlabel('Iterations');
      ylabel('Cost');

      subplot(2,4,[3;4;7;8])
      if i == 2, clf, end
      plot((0:size(H,2)-1)', - H(2*nthetag+2*p+2*1+2*q,:) + rbdo_parameters.target_beta, 'b.-' );
      title('\beta / Iterations')
      xlabel('Iterations');
      ylabel('\beta');

      savefigure = sprintf('hgsave N2LA.Temporary.Results.%d.fig;',i);
      eval(savefigure);

   end

   % Stopping criterias:

   if i > max_iter
      termflag = 0; % Max number of iterations overrun
   end

end

N2LAresults.number_of_iterations        = i-1;
N2LAresults.normalized_optimized_thetag = H(1:nthetag,:);
N2LAresults.optimized_thetag            = H(nthetag+1:2*nthetag,:);
N2LAresults.normalized_optimized_cost   = H(2*nthetag+p+1:2*nthetag+2*p,:)./very_initial_cost;
N2LAresults.optimized_cost              = H(2*nthetag+p+1:2*nthetag+2*p,:);
N2LAresults.optimized_beta              = - H(2*nthetag+2*p+2*1+2*q,:) + rbdo_parameters.target_beta;
N2LAresults.optimized_pf                = normcdf(-N2LAresults.optimized_beta);
N2LAresults.nfun                        = rbdo_nfun;
N2LAresults.nform                       = rbdo_nform;
N2LAresults.p                           = p;
N2LAresults.q                           = q;
N2LAresults.H                           = H;

if stepbystep_figures
   close(1);
end

figure

subplot(2,4,1)
hold on
plot((0:size(H,2)-1)', H(1:nthetag,:)','.-');
title('Normalized design variables (\theta_g) / Iterations')
xlabel('Iterations');
ylabel('Design variables (\theta_g)');
legend(char(gfundata(1).thetagname))

subplot(2,4,5)
hold on
plot((0:size(H,2)-1)', H(nthetag+1:2*nthetag,:)','.-');
title('Design variables (\theta_g) / Iterations')
xlabel('Iterations');
ylabel('Design variables (\theta_g)');
legend(char(gfundata(1).thetagname))

subplot(2,4,2)
plot((0:size(H,2)-1)', H(2*nthetag+p+1:2*nthetag+2*p,:)'./very_initial_cost, 'r.-' );
title('Normalized cost / Iterations')
xlabel('Iterations');
ylabel('Normalized cost');

subplot(2,4,6)
plot((0:size(H,2)-1)', H(2*nthetag+p+1:2*nthetag+2*p,:)', 'r.-' );
title('Cost / Iterations')
xlabel('Iterations');
ylabel('Cost');

subplot(2,4,[3;4;7;8])
plot((0:size(H,2)-1)', - H(2*nthetag+2*p+2*1+2*q,:) + rbdo_parameters.target_beta, 'b.-' );
title('\beta / Iterations')
xlabel('Iterations');
ylabel('\beta');