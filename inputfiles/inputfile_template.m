%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a template inputfile for FERUM 4.0.
% The data are given in at least 3 mandatory data structures, with various data fields.
%
% The 3 data structures are:
%  - probdata    ( stochastic model )
%  - gfundata    ( limit-state function )
%  - analysisopt ( analysis options )
%
% The number of data fields in these structures is variable and some of them may not have to be specified.
%
% If a global sensitivity analysis needs to be carried out (option 50 of FERUM 4.0), specific fields must be provided to analysisopt.
% See inputfile_Ellingwood_example44.m for details.
%
% If FERUM 4.0 makes use of an external code for evaluating the limit state function (e.g. a Finite Element code),
% the femodel data structure must be defined to set up environment variables and to pass arguments to the external code
%  - femodel     ( g-function based on an external code )
%
% If a RBDO analysis needs to be carried out (option 40 of FERUM 4.0), another data structure is required.
% It defines the RBDO parameters of the proposed N2LA analysis (cost, constraints, target reliability level, ...).
%  - rbdo_fundata ( cost and constraints )
%  - rbdo_parameters ( RBDO parameters )
%
% If a random field needs to be defined in the the stochastic model, the randomfield data structure is required.
% FERUM 4.0 does not provide yet any random field capabilities. See e.g. FERUMrandomfield of FERUM 3.1 version.
%  - randomfield ( random field model )
%
% The construction of this template file is based on the 4 following examples, presented in this order:
% - Example 1: inputfile_calrel_example4.m
% - Example 2: inputfile_nl_oscillators.m
% - Example 3: inputfile_SSLL116_geo.m ( limit-state function requiring calls to an external code, here code_Aster FE code )
% - Example 4: inputfile_rbdo_sphere_mmMPa.m ( RBDO analysis )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Clear possible old data in the Matlab workspace
clear probdata gfundata analysisopt femodel




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'PROBDATA'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.
% Use no spaces in names.
% probdata.name =  { 'name1' 'name2' ... } or { 'name1' 'name2' ... }'

%>>> Example 1 <<<<
probdata.name =  { 'X1'
                   'X2'
                   'X3'
                   'X4'
                   'X5'
                   'X6'};

%>>> Example 2 <<<<
probdata.name =  { 'mp'
                   'ms'
                   'kp'
                   'ks'
                   'zetap'
                   'zetas'
                   'Fs'
                   'S0' };

%>>> Example 3 <<<<
probdata.name =  { 'h'
                   'E'
                   'NU'
                   'F_VERT'
                   'EP1'
                   'RAY1' };

%>>> Example 4 <<<<
probdata.name =  { 'p0'
                   'r0'
                   'r1'
                   'fy' };
  
% Marginal distributions for each random variable
% probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];											  
%
% type: -1 = Parameter in reliability analysis (thetag)
%        0 = Deterministic parameter (cg)
%
%        1 = Normal distribution
%        2 = Lognormal distribution
%        3 = Gamma distribution
%        4 = Shifted exponential distribution
%        5 = Shifted Rayleigh distribution
%        6 = Uniform distribution
%        7 = Beta distribution
%        8 = Chi-square distribution
%
%       11 = Type I largest value distribution ( same as Gumbel distribution ) 
%       12 = Type I smallest value distribution
%       13 = Type II largest value distribution
%       14 = Type III smallest value distribution
%       15 = Gumbel distribution ( same as type I largest value distribution )  
%       16 = Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
%
%       18  (Reserved for Laplace distribution)
%       19  (Reserved for Pareto distribution)
%!      20  Generalized extremeve value distribution (GEV(p1 = shape, p2 = scale, p3 = location)) -> can be defined only with its parameters
%!      21  (Reserved for generalized Pareto distribution)
%!      25  Three-parameter lognormal distribution (p1 = shape, p2 = scale, p3 = thres)
%!      30  Custom probability distribution given with a sample -> kernel estimation is used in the analysis
%!      31  vector defined..
%
%       51 = Truncated normal marginal distribution
%
% Notes:
%
%    - Each field type, mean, stdv, startpoint, p1, p2, p3, p4, input_type must be filled in. If not used, input a dummy nan value.
%
%    - input_type = 0 when distributions are defined with mean and stdv only
%                   1 when distributions are defined with distribution parameters pi
%                    (i = 1, 2, 3 or 4, depending on the total number of parameters required)
%
%    - For the Type III smallest value marginal distribution, you must give the value of epsilon parameter as p3 when using the mean and stdv input (input_type = 0). 
%    - For the Beta marginal distribution , you must give the value of a parameter as p3 and b parameter as p4 when using the mean and stdv input (input_type = 0).
%    - For the Truncated normal marginal distribution, you must set input_type = 1 and give :
%         * the mean and stdv of the untruncated marginal distribution as p1 and p2 respectively
%         * the lower bound xmin and the upper bound xmax as p3 and p4 respectively
%
%    - startpoint stands for the starting point of the FORM analysis in the physical space
%
%!   - Sample based custom prob. distr. (30): 
%!        * p1 refers to 'tmp/custom_num2str(p1).mat' file where the sample is stored (under sample variable name!) 
%!        * p4 can be used to specify the power of powered custom distribution
%!        * other fields are irrelevant, if not given mean and std are calculated from inputs and startpoint is set to mean 
%!        * during the first call it generates the kernel distribution function: 'tmp/kernel_dist_num2str(p1).mat'
%!        * the sample can come from for example MCMC in case of Bayesian analysis
%
%    - Refer to ferum_pdf.m, ferum_cdf.m and distribution_parameter.m functions for more information on a specific distribution.

%>>> Example 1 <<<<
probdata.marg =  [  2    500   100   500  nan  nan  nan  nan 0 ;
                    2   2000   400  2000  nan  nan  nan  nan 0 ;
                    6      5   0.5     5  nan  nan  nan  nan 0 ;
                    2    450    90   450  nan  nan  nan  nan 0 ;
                    2   1800   360  1800  nan  nan  nan  nan 0 ;
                    6    4.5  0.45   4.5  nan  nan  nan  nan 0 ];

%>>> Example 2 <<<<
mu1 = 1.5;
Mu  = [ mu1 0.01 1.00 0.01 0.05 0.02 27.5 100 ]';
Cov = [ 0.1 0.1 0.2 0.2 0.4 0.5 0.1 0.1 ]';
probdata.marg =  [ 2*ones(8,1) Mu Mu.*Cov Mu  nan*zeros(8,4) zeros(8,1) ];

%>>> Example 3 <<<<
probdata.marg =    [ 0       0.1             0       0.1   nan   nan  nan  nan 0 ;
                     1  196200e6  196200e6*0.1  196200e6   nan   nan  nan  nan 0 ;
                     6       nan           nan       0.3  0.25  0.35  nan  nan 1 ;
                     1     -20e6     20e6*0.25     -20e6   nan   nan  nan  nan 0 ;
                     1      0.02      0.02*0.1      0.02   nan   nan  nan  nan 0 ;
                     1      0.05      0.05*0.1      0.05   nan   nan  nan  nan 0 ];

%>>> Example 4 <<<<
probdata.marg =  [  2   130     8   130 nan nan nan nan 0 ;
                   -1    50   nan    50 nan nan nan nan 0 ;
                   -1   100   nan   100 nan nan nan nan 0 ;
                    2   300    20   300 nan nan nan nan 0 ];



                  
% Correlation matrix
% This matrix is a square matrix with dimension equal to size(probdata.marg,1).
% Lines/columns corresponding to parameters in reliability analysis(thetag) or deterministic parameters (cg) are removed
% in a pre-processing stage of FERUM, by means of the update_data.m function.

%>>> Example 1 <<<<
probdata.correlation = [ 1.0 0.3 0.2 0.0 0.0 0.0 ;
                         0.3 1.0 0.2 0.0 0.0 0.0 ;
                         0.2 0.2 1.0 0.0 0.0 0.0 ;
                         0.0 0.0 0.0 1.0 0.3 0.2 ;
                         0.0 0.0 0.0 0.3 1.0 0.2 ;
                         0.0 0.0 0.0 0.2 0.2 1.0 ];

%>>> Example 2 <<<<
probdata.correlation = eye(8);

%>>> Example 3 <<<<
probdata.correlation = eye(6);

%>>> Example 4 <<<<
probdata.correlation = eye(4);


% Type of joint distribution
% transf_type = 1: jointly normal (no longer supported)
%               2: independent non-normal (no longer supported)
%               3: Nataf joint distribution (only available option)
%
probdata.transf_type = 3;

% Method for computation of the modified Nataf correlation matrix
% Ro_method = 0: use of approximations from ADK's paper (no longer supported)
%             1: exact, solved numerically
probdata.Ro_method   = 1;
                      
% Flag for computation of sensitivities w.r.t. means, standard deviations, parameters and correlation coefficients
% 1: all sensitivities assessed, 0: no sensitivities assessment
probdata.flag_sens   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'ANALYSISOPT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysisopt.echo_flag = 1;                   % 1: FERUM interactive mode, 0: FERUM silent mode

%>>> Example 1 <<<<

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
analysisopt.Recorded_u           = 1;        % 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
analysisopt.Recorded_x           = 1;        % 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations

% FORM, SORM analysis options
analysisopt.grad_flag            = 'ddm';    % 'ddm': direct differentiation, 'ffd': forward finite difference
analysisopt.ffdpara              = 1000;     % Parameter for computation of FFD estimates of gradients - Perturbation = stdv/analysisopt.ffdpara;
                                             % Recommended values: 1000 for basic limit-state functions, 50 for FE-based limit-state functions
analysisopt.ffdpara_thetag       = 1000;     % Parameter for computation of FFD estimates of dbeta_dthetag
                                             % perturbation = thetag/analysisopt.ffdpara_thetag if thetag ~= 0 or 1/analysisopt.ffdpara_thetag if thetag == 0;
                                             % Recommended values: 1000 for basic limit-state functions, 100 for FE-based limit-state functions

% Simulation analysis (MC,IS,DS,SS) and distribution analysis options
analysisopt.num_sim              = 1000;     % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)

% Simulation analysis (MC, IS) and distribution analysis options
analysisopt.sim_point            = 'origin'; % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.05;     % Target coefficient of variation for failure probability
analysisopt.lowRAM               = 0;        % 1: memory savings allowed, 0: no memory savings allowed

% Directional Simulation (DS) analysis options
analysisopt.dir_flag             = 'random'; % 'det': deterministic points uniformly distributed on the unit hypersphere using eq_point_set.m function
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


%>>> Example 2 <<<<

analysisopt.multi_proc           = 1;        % 1: block_size g-calls sent simultaneously
                                             %    - gfunbasic.m is used and a vectorized version of gfundata.expression is available.
                                             %      The number of g-calls sent simultaneously (block_size) depends on the memory
                                             %      available on the computer running FERUM.
                                             %    - gfunxxx.m user-specific g-function is used and able to handle block_size computations
                                             %      sent simultaneously, on a cluster of PCs or any other multiprocessor computer platform.
                                             % 0: g-calls sent sequentially
analysisopt.block_size           = 500000;   % Number of g-calls to be sent simultaneously

% FORM analysis options
analysisopt.i_max                = 500;      % Maximum number of iterations allowed in the search algorithm
analysisopt.e1                   = 0.001;    % Tolerance on how close design point is to limit-state surface
analysisopt.e2                   = 0.001;    % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code            = 0.025;    % 0: step size by Armijo rule, otherwise: given value is the step size
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
analysisopt.num_sim              = 10000;    % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)

% Simulation analysis (MC, IS) and distribution analysis options
analysisopt.sim_point            = 'origin'; % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.0001;   % Target coefficient of variation for failure probability
analysisopt.lowRAM               = 0;        % 1: memory savings allowed, 0: no memory savings allowed

% Directional Simulation (DS) analysis options
analysisopt.dir_flag             = 'random'; % 'det': deterministic points uniformly distributed on the unit hypersphere using eq_point_set.m function
                                             % 'random': random points uniformly distributed on the unit hypersphere
analysisopt.rho                  = 8;        % Max search radius in standard normal space for Directional Simulation analysis
analysisopt.tolx                 = 1e-5;     % Tolerance for searching zeros of g function
analysisopt.keep_a               = 0;        % Flag for storage of a-values which gives axes along which simulations are carried out
analysisopt.keep_r               = 0;        % Flag for storage of r-values for which g(r) = 0
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


%>>> Example 3 <<<<

analysisopt.multi_proc           = 0;        % 1: block_size g-calls sent simultaneously
                                             %    - gfunbasic.m is used and a vectorized version of gfundata.expression is available.
                                             %      The number of g-calls sent simultaneously (block_size) depends on the memory
                                             %      available on the computer running FERUM.
                                             %    - gfunxxx.m user-specific g-function is used and able to handle block_size computations
                                             %      sent simultaneously, on a cluster of PCs or any other multiprocessor computer platform.
                                             % 0: g-calls sent sequentially
analysisopt.block_size           = 1;        % Number of g-calls to be sent simultaneously
analysisopt.client_OS            = 'win';    % OS (Windows) running on the client

% FORM analysis options
analysisopt.i_max                = 10;       % Maximum number of iterations allowed in the search algorithm
analysisopt.e1                   = 0.01;     % Tolerance on how close design point is to limit-state surface
analysisopt.e2                   = 0.02;     % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code            = 0;        % 0: step size by Armijo rule, otherwise: given value is the step size
analysisopt.Recorded_u           = 1;        % 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
analysisopt.Recorded_x           = 1;        % 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations

% FORM, SORM analysis options
analysisopt.grad_flag            = 'ffd';    % 'ddm': direct differentiation, 'ffd': forward finite difference
analysisopt.ffdpara              = 50;       % Parameter for computation of FFD estimates of gradients - Perturbation = stdv/analysisopt.ffdpara;
                                             % Recommended values: 1000 for basic limit-state functions, 50 for FE-based limit-state functions
analysisopt.ffdpara_thetag       = 100;      % Parameter for computation of FFD estimates of dbeta_dthetag
                                             % perturbation = thetag/analysisopt.ffdpara_thetag if thetag ~= 0 or 1/analysisopt.ffdpara_thetag if thetag == 0;
                                             % Recommended values: 1000 for basic limit-state functions, 100 for FE-based limit-state functions

% Simulation analysis (MC,IS,DS,SS) and distribution analysis options
analysisopt.num_sim              = 200;      % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)

% Simulation analysis (MC, IS) and distribution analysis options
analysisopt.sim_point            = 'dspt';   % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.05;     % Target coefficient of variation for failure probability
analysisopt.lowRAM               = 0;        % 1: memory savings allowed, 0: no memory savings allowed

% Directional Simulation (DS) analysis options
analysisopt.dir_flag             = 'random'; % 'det': deterministic points uniformly distributed on the unit hypersphere using eq_point_set.m function
                                             % 'random': random points uniformly distributed on the unit hypersphere
analysisopt.rho                  = 7;        % Max search radius in standard normal space for Directional Simulation analysis
analysisopt.tolx                 = 5e-2;     % Tolerance for searching zeros of g function
analysisopt.keep_a               = 1;        % Flag for storage of a-values which gives axes along which simulations are carried out
analysisopt.keep_r               = 1;        % Flag for storage of r-values for which g(r) = 0
analysisopt.sigmafun_write2file  = 0;        % Set to 0

% Subset Simulation (SS) analysis options
analysisopt.width                = 2;        % Width of the proposal uniform pdfs
analysisopt.pf_target            = 0.1;      % Target probability for each subset step
analysisopt.flag_cov_pf_bounds   = 1;        % 1: calculate upper and lower bounds of the coefficient of variation of pf, 0: no calculation 
analysisopt.ss_restart_from_step = -1;       % i>=0 : restart from step i, -inf : all steps, no record (default), -1 : all steps, record all
analysisopt.flag_plot            = 0;        % 1: plots at each step (2 r.v. examples only), 0: no plots
analysisopt.flag_plot_gen        = 0;        % 1: intermediate plots for each MCMC chain (2 r.v. examples only), 0: no plots

% 2SMART analysis options
analysisopt.num_sim_sdu          = [ 50 50 ];         % Nu number of points ( [ first subset-like step, next steps], same values by default )
analysisopt.num_sim_norm         = [ 100 100 ];       % Nn number of points ( [ first subset-like step, next steps], same values by default )
nrv = size(probdata.marg,1);
analysisopt.buf_size             = round(1.25e7/nrv); % RAM-size dependent parameter (important for large number of random variables)
analysisopt.svm_buf_size         = 3500^2;            % RAM-size dependent parameter for eval_svm.m function
analysisopt.flag_var_rbf         = 0;                 % 0: cross validation at the 1st subset-like step only, 1: cross validation at all steps


%>>> Example 4 <<<<

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

%>>> Example 1 <<<<

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator  = 'basic';
gfundata(1).type       = 'expression';   % Do not change this field!

% Expression of the limit-state function:
gfundata(1).expression = 'c1 - X2./(1000*X3) - (X1./(200*X3)).^2 - X5./(1000*X6) - (X4./(200*X6)).^2';

% thetag parameters of the limit-state function (optional, may also be defined in probdata.marg)
% gfundata(1).thetag = [ (thetag1) (thetag2) ... ]
%                   or [ (thetag1) (thetag2) ... ]'
gfundata(1).thetag     = [ 1.7 ];

% Names of thetag parameters (thetag1, thetag2, ... if not defined):
% gfundata(1).thetagname = {'name1' 'name2' ... }
%                       or {'name1' 'name2' ... }'
gfundata(1).thetagname = { 'c1' };

% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens  = 1;


%>>> Example 2 <<<<

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator  = 'basic';
gfundata(1).type       = 'expression';   % Do no change this field!

% Expression of the limit-state function:
gfundata(1).expression = 'gfun_nl_oscillator(mp,ms,kp,ks,zetap,zetas,Fs,S0)';

% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens  = 0;


%>>> Example 3 <<<<

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator  = 'aster';
gfundata(1).type       = 'expression';   % Do no change this field!

% Expression of the limit-state function:
gfundata(1).expression = 'gext-u0';

% thetag parameters of the limit-state function (optional, may also be defined in probdata.marg)
% gfundata(1).thetag = [ (thetag1) (thetag2) ... ]
%                   or [ (thetag1) (thetag2) ... ]'
gfundata(1).thetag     = [];

% Names of thetag parameters (thetag1, thetag2, ... if not defined):
% gfundata(1).thetagname = {'name1' 'name2' ... }
%                       or {'name1' 'name2' ... }'
gfundata(1).thetagname = {};

% cg deterministic parameters of the limit-state function (optional, may also be defined in probdata.marg)
% gfundata(1).cg = [ (cg1) (cg2) ... ]
%               or [ (cg1) (cg2) ... ]'
gfundata(1).cg         = [ -0.6 ];

% Names of cg parameters (cg1, cg2, ... if not defined):
% gfundata(1).cgname = {'name1' 'name2' ... }
%                   or {'name1' 'name2' ... }'
gfundata(1).cgname     = { 'u0' };

% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens  = 1;


%>>> Example 4 <<<<


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

%>>> Example 3 <<<<

% Names of femodel.data variables (may be either realizations of random variables, deterministic or reliability parameters)
femodel.dataname = { 'h'
                     'E'
                     'NU'
                     'F_VERT'
                     'EP1'
                     'RAY1' };

% ID array to identify where the variables enter the FE model
% femodel.data = [ (type) (rank) (value) (parameter 1)(optional) (parameter 2)(optional)
%                  ... ];
%
% type:  -1 = parameter in reliability analysis (thetag)
%         0 = deterministic parameter (cg)
%         1 = Normal distribution
%         ... etc ... ( see field marg of probdata data structure )
%
% rank:  - rank in probdata.marg if the variable is defined in probdata.marg 
%        - nan if the variable is a deterministic or reliability parameter (cg or thetag respectively),
%          defined by gfundata data structure
%
% value: - nan if the variable is defined in probdata.marg
%          The value corresponds to the mean value of probdata.marg, if we are dealing with a deterministic 
%          or a reliability parameter (cg or thetag respectively), defined by probdata data structure (type = 0 or -1).
%          The value corresponds to a realization of the random variable, if type > 0.
%        - specify the value of the parameter if the corresponding rank is nan (usually, same value as
%          the one defined by gfundata data structure, i.e. gfundata.cg or gfundata.thetag)
%
% Remark:  femodel.data has priority over probdata.marg, gfundata(lsf).cg and gfundata(lsf).thetag.
%          Values of these data structures may be updated by those of femodel.data, if not identical,
%          e.g. type of probdata.marg, values of gfundata(lsf).cg or gfundata(lsf).thetag.
%          Before processing to a FE-based reliability analysis, make sure that the variables resulting
%          from the preprocessing stage of FERUM (result of the application of the update_data.m function
%          in ferum.m) are consistent with the analysis you want to carry out. 
%
femodel.data = [ 0    1  nan ;
                 1    2  nan ;
                 6    3  nan ;
                 1    4  nan ;
                 1    5  nan ;
                 1    6  nan ];

% Computer platform parameters - Enter paths with slashes (not backslashes), even under windows OS
femodel.mesh_generator      = 'geo';       % Use gmsh for mesh generation
femodel.mesh_generator_arg  = '-1';        % Argument passed to gmsh
femodel.jobdir_name         = 'exec';      % Computations are carried out in work_parent/execjob where job = 1, 2, ...
femodel.job_name            = 'ssll116b';  % Templates files processed and renamed as job_name.export, .geo, .comm
femodel.keep_optional_files = 1;           % 1: keep files in work_parent/execjob directories, 0: purge all directories
if strcmp(analysisopt.client_OS,'win')
   femodel.code_pathname    = 'C:/aster';                       % Code_Aster install directory
   femodel.script_path      = 'C:/temp/aster/seq1/scripts';     % Path of script files directory
                                                                % (content: call_aster.bat, call_gmsh.bat, extract_gext.awk)
   femodel.template_path    = 'C:/temp/aster/seq1/templates';   % Path of template files directory
                                                                % (content: template.geo, template.export, template.comm
   femodel.work_parent      = 'C:/temp/aster/seq1/execdir';     % Work parent directory
   femodel.client_dir       = 'C:/temp/aster/temp1';            % Temp directory
   femodel.gawk_path_name   = 'D:/gawk/gawk.exe';               % gawk.exe binary pathname (gawk used to postprocess result files)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'RBDO_FUNDATA' and 'RBDO_PARAMETERS'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%>>> Example 4 <<<<

                                                           % Cost function in the form: c_0 + c_f * p_f
rbdo_fundata.cost = {       '4/3.*pi.*(r1.^3 - r0.^3)'     % c_0 term
                      '0e2 * 4/3.*pi.*(r1.^3 - r0.^3)' };  % c_f term

rbdo_fundata.constraint = { '- r0 + 40'                    % Deterministic constraints: fi <= 0, i = 1,...,(q-1)
%                              'r1 - 150'
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