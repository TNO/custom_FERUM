function [ formresults, probdata ] = form_v3(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

% References
% NLopt     Steven G. Johnson, The NLopt nonlinear-optimization package, http://ab-initio.mit.edu/nlopt
%
%

% -------------------------------------------------------------------------
%  PREPROCESSING
% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
%  Constrained optimization
% -------------------------------------------------------------------------
% Bounds on x
bounds  = probdata.bounds;
if any(~isinf(bounds(:,1)))
    lb      = bounds(:,1);
else
    lb      = [];
end
if any(~isinf(bounds(:,2)))
    ub      = bounds(:,2);
else
    ub      = [];
end

% Starting point for the algorithm
x0      = marg(:,4);

form_solver = analysisopt.form_solver;

switch lower(form_solver)
    % Matlab Optimization Toolbox
    case {'interior-point', 'ip'}
        
        % The algorithm satisfies bounds at all iterations; Hessian based!
        alg                 = 'interior-point';
        options             = optimoptions('fmincon');
        options.Algorithm   = alg;
        options.ConstraintTolerance = 1e-2;
        options.OptimalityTolerance = 1e-2;
        options.Display     = 'iter';
        
        [xopt, fmin, exitflag, output] =  fmincon(@obj_fun,x0,[],[],[],[],lb,ub,@con_fun,options);
        
        % firstorderopt
%         Gopt                = output.firstorderopt;
        Gopt                = con_eq_fun(xopt);
%         nfun                = output.funcCount; % we might get rid of the
%         global nfun
        
        % Matlab Optimization Toolbox
    case 'sqp'
        
        alg                 = 'sqp'; % satisfies bounds at all iterations. The algorithm can recover from NaN or Inf results.
        options             = optimoptions('fmincon');
        options.Algorithm   = alg;
        options.ConstraintTolerance = 1e-2;
        options.OptimalityTolerance = 1e-2;
        options.Display     = 'iter';
        
        [xopt, fmin, exitflag, output] =  fmincon(@obj_fun,x0,[],[],[],[],lb,ub,@con_fun,options);
        
        Gopt                = output.constrviolation;
        nfun                = nfun + output.funcCount;
        
        % NLopt toolbox (free!), extended version to be able to handle bound
        % constrains in the form of linear constrains
    case 'cobyla'
%         x0 = x0';
        opt.algorithm       = NLOPT_LN_COBYLA;
        opt.min_objective   = @obj_fun;
        opt.h               = { (@(x) con_eq_fun(x))};
        opt.h_tol           = 1e-3;
        opt.xtol_rel        = 1e-4;
        opt.ftol_rel        = 1e-3;
%         opt.ftol_abs        = 1e-2;
        opt.maxeval         = 1e3;
        opt.verbose         = 1;
%         opt.initial_step    = 0.1*x0;
        
        % add bounds constrains if any random variables are bounded
%         ub = marg(:,2) + 7*marg(:,3);
        opt                 = cobyla_bound_con(lb, ub, opt);
        
        [xopt, fmin, exitflag] = nlopt_optimize(opt, x0);
        
        % wasteful but many results seem to be not accessible
        Gopt                = con_eq_fun(xopt);
        nfun                = nfun + NaN;

        % NLopt toolbox (free!)
    case 'slsqp'
        
        opt.algorithm       = NLOPT_LD_SLSQP;
        
    otherwise
        
        error(['Unknown FORM solver: ', form_solver])
end

% -------------------------------------------------------------------------
%  Post-processing
% -------------------------------------------------------------------------
uopt                    = x_to_u(xopt,probdata);

formresults.beta        = fmin;
formresults.dsptx       = xopt;
formresults.dsptu       = uopt;
formresults.G           = Gopt;
formresults.nfun        = nfun;
formresults.solver      = form_solver;
formresults.exitflag    = exitflag;

% -------------------------------------------------------------------------
% NESTED FUNCTIONS - for the constrained optmization
% -------------------------------------------------------------------------
    function b = obj_fun(x)
        x = x(:);
        u = x_to_u(x,probdata);
        b = norm(u);
    end

    function [c, G] = con_fun(x)
        c = con_ineq_fun(x);
        G = con_eq_fun(x);
    end

    function G = con_eq_fun(x)
        x = x(:);
        G = gfun(lsf,x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
    end

% dummy - needed because of how Matlab's fmincon works
    function c = con_ineq_fun(~)
        c = -1;
    end

end



