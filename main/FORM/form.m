function [ formresults, probdata ] = form(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

marg        = probdata.marg;
R           = probdata.correlation;
transf_type = probdata.transf_type;
flag_sens   = probdata.flag_sens;

if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end

i_max     = analysisopt.i_max;
e1        = analysisopt.e1;
e2        = analysisopt.e2;
step_code = analysisopt.step_code;
grad_flag = lower(analysisopt.grad_flag);
if isfield(analysisopt,'Recorded_u')
   if analysisopt.Recorded_u
      Recorded_u_flag = 1;
      Recorded_u = [];
   else
      Recorded_u_flag = 0;
   end
else
   Recorded_u_flag = 0;    
end
if isfield(analysisopt,'Recorded_x')
   if analysisopt.Recorded_x
      Recorded_x_flag = 1;
      Recorded_x = [];
   else
      Recorded_x_flag = 0;    
   end
else
   Recorded_x_flag = 0;    
end
if isfield(analysisopt,'Recorded_sens')
   if analysisopt.Recorded_sens
      Recorded_sens_flag = 1;
      Recorded_grad_G = [];
      Recorded_alpha  = [];
      Recorded_imptg  = [];
   else
      Recorded_sens_flag = 0;    
   end
else
   Recorded_sens_flag = 0;    
end

% Modify correlation matrix and perform Cholesky decomposition
if ~isfield(probdata,'Lo')
    
   if transf_type == 3
    
      % Compute corrected correlation coefficients
      switch probdata.Ro_method
         case 1
            if echo_flag
               disp([' '])
               disp('Computation of modified correlation matrix R0')
               disp('Takes some time if sensitivities are to be computed with gamma (3), beta (7) or chi-square (8) distributions.')
               disp('Please wait... (Ctrl+C breaks)')
            end
            [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , flag_sens);
      end     
      probdata.Ro = Ro;
   
      % Cholesky decomposition
      [Lo,ierr] = my_chol(Ro); if  ierr>0, return, end
      probdata.Lo = Lo;
   
      iLo = inv(Lo);
      probdata.iLo = iLo;

   end
   
else
    
   Ro = probdata.Ro;
   if isfield(probdata,'dRo_drho')
      dRo_drho = probdata.dRo_drho;
      dRo_dthetafi = probdata.dRo_dthetafi;
      dRo_dthetafj = probdata.dRo_dthetafj;
   elseif flag_sens
      % Compute corrected correlation coefficients
      switch probdata.Ro_method
         case 1
            if echo_flag
               disp([' '])
               disp('Computation of modified correlation matrix R0')
               disp('Takes some time if sensitivities are to be computed with gamma (3), beta (7) or chi-square (8) distributions.')
               disp('Please wait... (Ctrl+C breaks)')
            end
            [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , flag_sens);
      end     
   end

   Lo = probdata.Lo;
   iLo = probdata.iLo;

end

% Compute starting point for the algorithm
x = marg(:,4);
u = x_to_u(x,probdata);

if any(isnan(x))
   keyboard
end

% Set parameters for the iterative loop
i = 1;         % Initialize counter
conv_flag = 0; % Convergence is achieved when this flag is set to 1
shit = 0;
% Perform iterative loop to find the design point
while conv_flag == 0;
%     keyboard
   if echo_flag
      disp('..................................')
      disp('Now carrying out iteration number:'),disp(i)
   end

   % Transformation from u to x space
   x = u_to_x(u,probdata);
%    if any(isnan(x))
%        keyboard
%    end
   J_u_x = jacobian(x,u,probdata);
   J_x_u = inv(J_u_x);
   
   %=======================================================================
   % Stop FORM iteration and return NaN beta if ill-conditioned (close to singular) matrix is detected
   id = warning('query','last');
   if ~isempty(id)
       last_warn = id.identifier;
       if any(strcmp(last_warn, 'MATLAB:illConditionedMatrix')) && i > 1
           conv_flag = 1;
           % overwrite last warning and provide info about what is happening
           warning('MyCalibration:illconditionedFORM', 'Aborting FORM iteration and returning beta = NaN (hopefully :)). No convergence!')
           shit = 1;
       end
   end
   %=======================================================================
   
   % Evaluate limit-state function and its gradient
   [ G, grad_g ] = gfun(lsf,x,grad_flag,probdata,analysisopt,gfundata,femodel,randomfield);
   
   grad_G = (grad_g' * J_x_u)';
   Recorded_G_function_values(i) = G;
   
   % Set scale parameter Go and inform about struct. resp.
   if i == 1
      Go = G;
      if echo_flag
         disp('Value of limit-state function in the first step:')
         disp(G)
      end
   end
   
   % Compute alpha vector
   alpha = -grad_G / norm(grad_G);                                      % Alpha vector
   % Compute gamma vector
   D_prime = diag(diag(   sqrt( J_x_u * J_x_u' )    ));                 % (intermediate computation)
    imptg = (alpha'*J_u_x*D_prime/norm(alpha'*J_u_x*D_prime))';          % Importance vector gamma

   % Check convergence
   if ( (abs(G/Go)<e1) & (norm(u-alpha'*u*alpha)<e2) ) | i==i_max
      conv_flag = 1;
   end
   Recorded_G_criterion(i) = abs(G/Go);
   Recorded_u_criterion(i) = norm(u-alpha'*u*alpha);
   if Recorded_u_flag
      Recorded_u = [ Recorded_u u ];
   end
   if Recorded_x_flag
      Recorded_x = [ Recorded_x x ];
   end
	
   if Recorded_sens_flag
      Recorded_grad_G = [ Recorded_grad_G grad_G ];                      % Value of grad_G(u)
      Recorded_alpha  = [ Recorded_alpha alpha ];
      Recorded_imptg  = [ Recorded_imptg imptg ];
 
   end

   % Take a step if convergence is not achieved
   if conv_flag == 0;
      
      % Determine search direction
      d = search_dir(G,grad_G,u);
%       if any(isnan(d))
%           keyboard
%       end
      % Determine step size
      if step_code == 0
% %           disp('step')
         step = step_size_multiproc(lsf,G,grad_G,u,d,probdata,analysisopt,gfundata,femodel,randomfield);
      else
         step = step_code;
      end
      Recorded_step_size_values(i) = step;
      
      %!!!!!!!!!!!!!!! WARNING!
      if norm(d) > 1
        d = d./norm(d);
      end
      %!!!!!!!!!!!!!!! WARNING!
      
      % Determine new trial point
      u_new = u + step * d;
      
      
      % Prepare for a new round in the loop
      u = u_new;
      i = i + 1;
      
      Recorded_beta(i) = alpha' * u;
% %       disp(['beta=',num2str(alpha' * u)])
   end
   
end   

beta = alpha' * u;


if echo_flag & Recorded_u_flag
   nrv = size(marg,1);
   if nrv == 2
      plot(Recorded_u(1,:),Recorded_u(2,:),'r*-')
   end
   if nrv == 3
      plot3(Recorded_u(1,:),Recorded_u(2,:),Recorded_u(3,:),'r*-')
   end
end


%  Post-processing 
formresults.iter   = i;                                                 % Number_of_iterations
beta               = alpha' * u;
formresults.beta   = beta;                                              % Reliability index beta
formresults.pf1    = normcdf(-formresults.beta);                        % Probability of failure pf1
% formresults.pf1    = normcdf_mp(-formresults.beta);                        % Probability of failure pf1
formresults.dsptu  = u;                                                 % Design point u_star
formresults.G      = G;                                                 % Value of G(u_star)
formresults.grad_G = grad_G;                                            % Value of grad_G(u_star)
formresults.alpha  = alpha;                                             % Alpha vector
formresults.dsptx  = x;                                                 % Design point in original space
formresults.imptg  = imptg;                                             % Importance vector gamma
formresults.gfcn   = Recorded_G_function_values;                        % Recorded values of the limit-state function
formresults.stpsz  = Recorded_step_size_values;                         % Recorded step size values
formresults.e1     = Recorded_G_criterion;                              % Recorded values of g-criterion
formresults.e2     = Recorded_u_criterion;
formresults.Recorded_beta     = Recorded_beta;
% Recorded values og u-criterion
if Recorded_u_flag
   formresults.Recorded_u = Recorded_u;                                 % Recorded values of u
end
if Recorded_x_flag
   formresults.Recorded_x = Recorded_x;                                 % Recorded values of x
end
if Recorded_sens_flag
   formresults.Recorded_grad_G = Recorded_grad_G;
   formresults.Recorded_alpha  = Recorded_alpha;
   formresults.Recorded_imptg  = Recorded_imptg;
end

% shit = 1; % always show the plot
if i == i_max
    warning('Maximum number of iterations was reached before convergence.')
end

% if i == i_max || shit == 1  
%    if analysisopt.diagn_plot == 1
%        figure('Position', [100, -200, 800, 800])
%        max_istep = [analysisopt.i_max, 20];
%        for jj = 1:2
%            subplot(3,2,1+1*(jj-1))
%            ll = min(i-1, max_istep(jj));
%            tmp = Recorded_x(:,end-ll:end);
%            tt = bsxfun(@rdivide, tmp, tmp(:,1));
%            plot(tt')
%            ylabel(['normalized Recorded x (last ', num2str(ll),')'])
% 
%            subplot(3,2,3+1*(jj-1))
%            plot(Recorded_beta(end-ll:end))
%            ylabel(['Recorded beta (last ', num2str(ll),')'])
% 
%            subplot(3,2,5+1*(jj-1))
%            plot(Recorded_G_criterion(end-ll:end))
%            hold on
%            plot(Recorded_u_criterion(end-ll:end))
%            ylabel(['Recorded criteria (last ', num2str(ll),')'])
%        end
% %    disp('Press a key to continue!')
% %    pause;
%    end
%    
% end


if probdata.flag_sens
   switch transf_type
      case 3
         [ dbeta_drho, dpf_drho, dbeta_dthetaf, dpf_dthetaf, dbeta_dthetaf_approx, dbeta_dthetaf_ratio ] = ...
             sensitivities_wrt_thetaf( beta, alpha, x, u, marg, Ro, Lo, iLo, dRo_drho, dRo_dthetafi, dRo_dthetafj );
         delta = dbeta_dthetaf(:,1) .* marg(:,3);
         eta = dbeta_dthetaf(:,2) .* marg(:,3);
         formresults.dbeta_dthetaf        = dbeta_dthetaf;              % Sensitivities of beta index w.r.t. to distribution parameters
         formresults.dbeta_dthetaf_approx = dbeta_dthetaf_approx;       % Approximated sensitivities of beta index w.r.t. to distribution parameters (usual approximation)
         formresults.dbeta_dthetaf_ratio  = dbeta_dthetaf_ratio;        % Ratio ( alpha' * diLo_dthetaf * z ) / ( alpha' * iLo * J_z_thetaf(:,6*(i-1)+k) ), k = 1 to 6
         formresults.dpf_dthetaf          = dpf_dthetaf;                % Sensitivities of failure probability w.r.t. to distribution parameters
         formresults.dbeta_drho           = dbeta_drho;                 % Sensitivities of beta index w.r.t. to correlation coefficients
         formresults.dpf_drho             = dpf_drho;                   % Sensitivities of failure probability w.r.t. to correlation coefficients
         formresults.delta                = delta;                      % "Normalized" sensitivities of beta index w.r.t. mean
         formresults.eta                  = eta;                        % "Normalized" sensitivities of beta index w.r.t. standard deviation
   end
end


if gfundata(lsf).flag_sens
   if isfield(gfundata(1),'thetag')
      [ dbeta_dthetag, dpf_dthetag ] = sensitivities_wrt_thetag(lsf,beta,x,G,grad_G,probdata,analysisopt,gfundata,femodel,randomfield);
   else
      dbeta_dthetag = [];
      dpf_dthetag = [];
   end
   formresults.dbeta_dthetag = dbeta_dthetag;                           % Sensisitivities of beta index w.r.t. limit-state parameters
   formresults.dpf_dthetag = dpf_dthetag;                               % Sensisitivities of failure probability w.r.t. limit-state parameters
end
   
   
formresults.nfun = nfun;