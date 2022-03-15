function [ ALLformresults, probdata, analysisopt, gfundata ] = form_multiple_dspts(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

% Set method specific constants
epsil = 0.5;  % Recommended value: 0.2-0.5 
gamm  = 1.1;
delt  = 0.75;

% Forces gradients to be computed by finite differences
% DDM version would need the differentiation of Bi's w.r.t. x 
analysisopt.grad_flag = 'ffd';

% Find number of random variables
nrv = size(probdata.marg,1);

marg        = probdata.marg;
R           = probdata.correlation;
transf_type = probdata.transf_type;
flag_sens   = probdata.flag_sens;

if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end

% Modify correlation matrix and perform Cholesky decomposition
if ~isfield(probdata,'Lo')
    
   if transf_type == 3
    
      % Compute corrected correlation coefficients
      switch probdata.Ro_method
         case 0
            Ro = mod_corr( marg, R );
         case 1
            disp([' '])
            disp('Computation of modified correlation matrix R0')
            disp('Takes some time if sensitivities are to be computed with gamma (3), beta (7) or chi-square (8) distributions.')
            disp('Please wait... (Ctrl+C breaks)')
            [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , flag_sens);
      end     
      probdata.Ro = Ro;
      probdata.dRo_drho = dRo_drho
      probdata.dRo_dthetafi = dRo_dthetafi;
      probdata.dRo_dthetafj = dRo_dthetafj;
   
      % Cholesky decomposition
      [Lo,ierr] = my_chol(Ro); if  ierr>0, return, end
      probdata.Lo = Lo;
   
      iLo = inv(Lo);
      probdata.iLo = iLo;

   end

end

% Initializations
i     = 0;
idspt = 0;

while 1

   % Set starting point
   if i > 0
      u0 = -epsil * sum( dsptu, 2 );
      x0 = u_to_x(u0,probdata);   
      probdata.marg(:,4) = x0;
   end
      
   [ formresults, probdata ] = form(lsf,probdata,analysisopt,gfundata,femodel,randomfield);

   if i > 0
      
      % Test if solution is not a spurious one
      dist_found_to_previous = dot( formresults.dsptu*ones(1,idspt) - dsptu, formresults.dsptu*ones(1,idspt) - dsptu ).^0.5;
      i_in_bulge = find( dist_found_to_previous < r );
      
      if isempty(i_in_bulge) & ( formresults.iter < analysisopt.i_max )
         
         idspt = idspt + 1;
         pf1(idspt) = formresults.pf1;
         
         if echo_flag
            fprintf(1,'\nDesign point #%d\n', idspt);
            fprintf('\nbeta = %.4f\n', formresults.beta);
            fprintf('\npf   = %.3e\n', pf1(idspt));
            fprintf('\nu*   = \n');
            fprintf('%13.4f\n', formresults.dsptu);
            fprintf(1,'\npf/sum(pf) = %.4f\n', pf1/sum(pf1));
            answer = input('\nNew search ? y / n\n', 's');
            while ~( strcmp( lower(answer) ,'y') | strcmp( lower(answer) ,'n') )
               fprintf(1,'Incorrect answer, type y or n\n');
               answer = input('New search ? y / n\n', 's');
            end
            fprintf('\n\n\n');
         else
            answer = 'y';
         end
         
         r(idspt) = gamm * formresults.beta;
         dsptu(1:nrv,idspt) = formresults.dsptu;
         ALLformresults(idspt) = formresults;
         
         if strcmp( lower(answer) ,'n')
            break
         end
         
      elseif ~isempty(i_in_bulge)
         
         if echo_flag
            fprintf(1,'\nPotential design point #%d falls in bulge of design point #%d\n', [ (idspt+1)*ones(1,length(i_in_bulge)) i_in_bulge])
            fprintf(1,'\nSpurious design point #%d\n', idspt+1);
            fprintf('\nbeta = %.4f\n', formresults.beta);
            fprintf('\nu* = \n');
            fprintf('%13.4f\n', formresults.dsptu);
            answer = input('\nNew search ? y / n\n', 's');
            while ~( strcmp( lower(answer) ,'y') | strcmp( lower(answer) ,'n') )
               fprintf(1,'Incorrect answer, type y or n\n');
               answer = input('New search ? y / n\n', 's');
            end
         else
            answer = 'n';
         end
         
         if strcmp( lower(answer) ,'n')
            break
         end
         
      elseif formresults.iter == analysisopt.i_max
         
         if echo_flag
            fprintf(1,'\nMax number of iterations reached in iHL-RF algoritm\n');
            fprintf(1,'\nSpurious design point #%d\n', idspt+1);
            fprintf('\nbeta = %.4f\n', formresults.beta);
            fprintf('\nu* = \n');
            fprintf('%13.4f\n', formresults.dsptu);
            answer = input('\nNew search ? y / n\n', 's');
            while ~( strcmp( lower(answer) ,'y') | strcmp( lower(answer) ,'n') )
               fprintf(1,'Incorrect answer, type y or n\n');
               answer = input('New search ? y / n\n', 's');
            end
         else
            answer = 'n';
         end
         
         if strcmp( lower(answer) ,'n')
            break
         end
         
      end
      
   else
      
      idspt = idspt + 1;
      pf1(idspt) = formresults.pf1;

      if echo_flag
         fprintf(1,'\nDesign point #%d\n', idspt);
         fprintf(1,'\nbeta = %.4f\n', formresults.beta);
         fprintf(1,'\npf   = %.3e\n', pf1(idspt));
         fprintf(1,'\nu*   = \n');
         fprintf(1,'%13.4f\n', formresults.dsptu);
         fprintf(1,'\npf/sum(pf) = %.4f\n', pf1/sum(pf1));
         fprintf('\nPress a key to continue');
         pause
         fprintf('\n\n\n');
      end
      
      r(idspt) = gamm * formresults.beta;
      dsptu(1:nrv,idspt) = formresults.dsptu;
      ALLformresults(idspt) = formresults;
      
   end
   
   i = i + 1;
   
   betai(i) = formresults.beta;
   dpstui(1:nrv,i) = formresults.dsptu;
   norm_of_grad_Gi(i) = norm(formresults.grad_G);
   si(i) = delt * betai(i) * norm_of_grad_Gi(i) / ( ( gamm * betai(i) )^2 - ( delt * betai(i) )^2 )^2;
   ri(i) = gamm * betai(i);
   
   gfundata(lsf).bulge = 1;

   gfundata(lsf).m_minus_1 = i;
   gfundata(lsf).si(i) = si(i);
   gfundata(lsf).ri(i) = ri(i);
   gfundata(lsf).dsptui(1:nrv,i) = dpstui(1:nrv,i);
   
   eval(['save gfundata_' num2str(i) '.mat gfundata']);
   
end

% Restore original g function (but keep Bi's information)
gfundata(lsf).bulge = 0;

% Store ALLformresults in analysisopt structure
analysisopt.ALLformresults = ALLformresults;