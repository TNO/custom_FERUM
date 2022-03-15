%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'RBDO_FUNDATA' and 'RBDO_PARAMETERS'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                           % Cost function in the form: c_0 + c_f * p_f
rbdo_fundata.cost = {       '4/3.*pi.*(r1.^3 - r0.^3)'     % c_0 term
                      '0e2 * 4/3.*pi.*(r1.^3 - r0.^3)' };  % c_f term

rbdo_fundata.constraint = { '- r0 + 40e-3'                    % Deterministic constraints: fi <= 0, i = 1,...,(q-1)
%                              'r1 - 150e-3'
                              'r0 - r1'
                          };

rbdo_parameters.alpha                    = 0.5;            % ph_quadprog parameter
rbdo_parameters.beta                     = 0.6;            % ph_quadprog parameter
rbdo_parameters.gamma                    = 3;              % ph_quadprog parameter
rbdo_parameters.delta                    = 1;              % ph_quadprog parameter

rbdo_parameters.steplim                  = 50;             % Max number of steps in stepsize calculation (see ph_quadprog.m)

rbdo_parameters.max_iter                 = 7;              % Max number of iterations of N2LA algorithm
rbdo_parameters.target_beta              = 5;              % Target beta reliability index
rbdo_parameters.method                   = 'FORM';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'PROBDATA'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.
probdata.name =  { 'p0'
                   'r0'
                   'r1'
                   'fy' };

% Marginal distributions for each random variable
% probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];
probdata.marg =  [  2   130e6   8e6   130e6 nan nan nan nan 0 ;
                   -1   50e-3   nan   50e-3 nan nan nan nan 0 ;
                   -1  100e-3   nan  100e-3 nan nan nan nan 0 ;
                    2   300e6  20e6   300e6 nan nan nan nan 0 ];

% Correlation matrix
probdata.correlation = eye(4);

probdata.transf_type = 3;
probdata.Ro_method   = 1;
probdata.flag_sens   = 1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'ANALYSISOPT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysisopt.multi_proc           = 1;        % 1: block_size g-calls sent simultaneously
                                             %    - gfunbasic.m is used and a vectorized version of gfundata.expression is available.
                                             %      The number of g-calls sent simultaneously (block_size) depends on the memory
                                             %      available on the computer running FERUM.
                                             %    - gfunxxx.m user-specific g-function is used and able to handle block_size computations
                                             %      sent simultaneously, on a cluster of PCs or any other multiprocessor computer platform.
                                             % 0: g-calls sent sequentially
analysisopt.block_size           = 500000;   % Number of g-calls to be sent simultaneously

% FORM analysis options
analysisopt.i_max                = 100;      % Maximum number of iterations allowed in the search algorithm
analysisopt.e1                   = 0.001;    % Tolerance on how close design point is to limit-state surface
analysisopt.e2                   = 0.001;    % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code            = 0;        % 0: step size by Armijo rule, otherwise: given value is the step size
analysisopt.Recorded_u           = 0;        % 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
analysisopt.Recorded_x           = 0;        % 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations

% FORM, SORM analysis options
analysisopt.grad_flag            = 'ffd';    % 'ddm': direct differentiation, 'ffd': forward finite difference
analysisopt.ffdpara              = 1000;     % Parameter for computation of FFD estimates of gradients - Perturbation = stdv/analysisopt.ffdpara;
                                             % Recommended values: 1000 for basic limit-state functions, 50 for FE-based limit-state functions
analysisopt.ffdpara_thetag       = 1000;     % Parameter for computation of FFD estimates of dbeta_dthetag
                                             % perturbation = thetag/analysisopt.ffdpara_thetag if thetag ~= 0 or 1/analysisopt.ffdpara_thetag if thetag == 0;
                                             % Recommended values: 1000 for basic limit-state functions, 100 for FE-based limit-state functions

% Simulation analysis (MC,IS,DS,SS) and distribution analysis options
analysisopt.num_sim              = 1e4;      % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)

% Simulation analysis (MC, IS) and distribution analysis options
analysisopt.sim_point            = 'origin'; % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.05;     % Target coefficient of variation for failure probability
analysisopt.lowRAM               = 0;        % 1: memory savings allowed, 0: no memory savings allowed

% Directional Simulation (DS) analysis options
analysisopt.dir_flag             = 'det';    % 'det': deterministic points uniformly distributed on the unit hypersphere using eq_point_set.m function
                                             % 'random': random points uniformly distributed on the unit hypersphere
analysisopt.rho                  = 8;        % Max search radius in standard normal space for Directional Simulation analysis
analysisopt.tolx                 = 1e-5;     % Tolerance for searching zeros of g function
analysisopt.keep_a               = 1;        % Flag for storage of a-values which gives axes along which simulations are carried out
analysisopt.keep_r               = 1;        % Flag for storage of r-values for which g(r) = 0
analysisopt.sigmafun_write2file  = 0;        % Set to 0

% Subset Simulation (SS) analysis options
analysisopt.width                = 2;        % Width of the proposal uniform pdfs
analysisopt.pf_target            = 0.1;      % Target probability for each subset step
analysisopt.flag_cov_pf_bounds   = 1;        % 1: calculate upper and lower bounds of the coefficient of variation of pf, 0: no calculation 
analysisopt.ss_restart_from_step = -inf;     % i>=0 : restart from step i, -inf : all steps, no record (default), -1 : all steps, record all
analysisopt.flag_plot            = 0;        % 1: plots at each step (2 r.v. examples only), 0: no plots
analysisopt.flag_plot_gen        = 0;        % 1: intermediate plots for each MCMC chain (2 r.v. examples only), 0: no plots

% 2SMART analysis options
analysisopt.num_sim_sdu          = [ 50 50 ];         % Nu number of points ( [ first subset-like step, next steps], same values by default )
analysisopt.num_sim_norm         = [ 100 100 ];       % Nn number of points ( [ first subset-like step, next steps], same values by default )
nrv = size(probdata.marg,1);
analysisopt.buf_size             = round(1.25e7/nrv); % RAM-size dependent parameter (important for large number of random variables)
analysisopt.svm_buf_size         = 3500^2;            % RAM-size dependent parameter for eval_svm.m function
analysisopt.flag_var_rbf         = 0;                 % 0: cross validation at the 1st subset-like step only, 1: cross validation at all steps




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'GFUNDATA' (one structure per gfun)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator  = 'basic';
gfundata(1).type       = 'expression';   % Do no change this field!

% Expression of the limit-state function:
gfundata(1).expression = 'fy - 1.5*p0.*(r1.^3)./((r1.^3)-(r0.^3))';

% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens  = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'FEMODEL'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

femodel = [];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'RANDOMFIELD'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomfield = [];