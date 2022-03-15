function [ lambdaresults, probdata ] = lambda_method(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

% TODO: 
% 1)using the FERUM infrastructure to evaluate the gradient unneccesary
% limit state function evaluations are performed for 'ffd' it is small
% for 'ddm' considerable
% 2) still the form method seems to do some strange gfuneval counting -
% should be checked
%   simple example for benchmarking: D:\Working folder\Matlab working folder\lambda_method\myexample
%
%

% warning('For lambda method the limit state function should have the form of g(e,x) = e - R(x)')
global nfun
nfun = 0;

%Parameters of the lambda method
%!!!!!!!!!!!!!!
lambda  = analysisopt.lambda; 
m       = 1; % exponent
%!!!!!!!!!!!!!!

marg        = probdata.marg;
R           = probdata.correlation;
transf_type = probdata.transf_type;
flag_sens   = probdata.flag_sens;

if isfield(analysisopt,'echo_flag')
    echo_flag = analysisopt.echo_flag;
else
    echo_flag = 1;
end

grad_flag = lower(analysisopt.grad_flag);

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
xx = marg(:,4);
ini = x_to_u(xx, probdata);

% keyboard
% fulfill the collinearity condition: system of nonlinear equations, Broyden's method with Armijo rule 
[uu, ~, ~, iter] = brsola(ini, @collinear_crit, [1e-4, 1e-4]); % WARNING!

% Compute alpha vector
alpha = -grad_G / norm(grad_G);

beta = alpha' * uu;

r = gfun(lsf,xx,grad_flag,probdata,analysisopt,gfundata,femodel,randomfield);

% Nested functions for lambda-method
    function grad_uR = grad_uR_fun(uu)
        % Transformation from u to x space
        xx = u_to_x(uu,probdata);
        
        J_u_x = jacobian(xx,uu,probdata);
        J_x_u = inv(J_u_x);
        
        % Evaluate limit-state function and its gradient
        [ ~, grad_g ] = gfun(lsf,xx,grad_flag,probdata,analysisopt,gfundata,femodel,randomfield);
        
        grad_G  = (grad_g' * J_x_u)';
        grad_uR = -grad_G;
    end

    function F = collinear_crit(uu)
        grad_uR = grad_uR_fun(uu);
        F =  uu - lambda/(norm(grad_uR)^m)*grad_uR;
    end


%  Post-processing
lambdaresults.r      = r;                                                 % Threshold value
lambdaresults.iter   = iter;                                              % Number of iterations
lambdaresults.beta   = beta;                                              % Reliability index beta
lambdaresults.pf1    = normcdf(-beta);                                    % Probability of failure pf1
lambdaresults.dsptu  = uu;                                                % Design point u_star
lambdaresults.grad_G = grad_G;                                            % Value of grad_G(u_star)
lambdaresults.dsptx  = xx;                                                % Design point in physical space
lambdaresults.alpha  = alpha;                                             % Sensitivity vector
lambdaresults.nfun   = nfun;                                              % Number of gfun evaluations

end