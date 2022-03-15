% rigid_frame_system.m reliability assessment
%
% [Zhong, Zhao, Ang (2007). A method for computing reliabiliy bound of series structural systems.
%  Reliability and Optimization of Structural Systems: Assessment, Desing and Life-Cycle Perfornance. London]
%  web: http://books.google.hu/books?id=2OCWkLicuuEC&pg=PA259&lpg=PA259&dq=A+method+for+computing+reliability+bound+of+series+structural+systems&source=bl&ots=FJQ_EY2GJE&sig=jIe3NvIAvq_J3dO7FhMJYAjUwHk&hl=en&sa=X&ei=yv71U-6qBOuy7AaOn4BI&redir_esc=y#v=onepage&q=A%20method%20for%20computing%20reliability%20bound%20of%20series%20structural%20systems&f=false
% [crossref: Ono, T., Idota, H., and Dozuka, A. (1990). “Reliability evaluation of structural systems using higher-order
%  moment standardization technique.” J. of Struct. Constr. Engng., AIJ, No. 418, 71–79 (In Japanese).]
%
% results from the paper:
%              1      2      3      4      5      6
% beta_comp = [3.247, 3.551, 3.562, 3.562, 3.784, 3.848]
%
%!!??[mm], [N], [s], [kg]??!!



clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DATA FIELDS IN 'PROBDATA'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PROBABILISITC MODEL - RANDOM VARIABLES
% resistances
M1_m                        = 500;              % mean value (only mean or char should be defined!)
M1_cov                      = 75/500;           % coefficient of variation
M1_dist                     = 2;                % distribution type
M1_input                    = 0;                % input type
M1_std                      = M1_m*M1_cov;

M2_m                        = 500;              % mean value (only mean or char should be defined!)
M2_cov                      = 75/500;           % coefficient of variation
M2_dist                     = 2;                % distribution type
M2_input                    = 0;                % input type
M2_std                      = M2_m*M2_cov;

M3_m                        = 667;              % mean value (only mean or char should be defined!)
M3_cov                      = 100/667;          % coefficient of variation
M3_dist                     = 2;                % distribution type
M3_input                    = 0;                % input type
M3_std                      = M3_m*M3_cov;

% effect/action
S1_m                         = 50;               % mean value (only mean or char should be defined!)
S1_cov                       = 15/50;            % coefficient of variation
S1_dist                      = 2;                % distribution type
S1_input                     = 0;                % input type
S1_std                       = S1_m*S1_cov;

S2_m                         = 100;              % mean value (only mean or char should be defined!)
S2_cov                       = 10/100;           % coefficient of variation
S2_dist                      = 2;                % distribution type
S2_input                     = 0;                % input type
S2_std                       = S2_m*S2_cov;


% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.
probdata.name = { 'M1'
    'M2'
    'M3'
    'S1'
    'S2'
    'als'};

% Marginal distributions for each random variable
% probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];
probdata.marg =  [M1_dist, M1_m, M1_std, M1_m,  NaN, NaN, NaN, NaN, 0;
    M2_dist, M2_m, M2_std, M2_m,  NaN, NaN, NaN, NaN, 0;
    M3_dist, M3_m, M3_std, M3_m,  NaN, NaN, NaN, NaN, 0;
    S1_dist, S1_m, S1_std, S1_m,  NaN, NaN, NaN, NaN, 0;
    S2_dist, S2_m, S2_std, S2_m,  NaN, NaN, NaN, NaN, 0;
    0,       1,    0,      NaN,   NaN, NaN, NaN, NaN, 0];

% Correlation matrix
probdata.correlation = eye(length(probdata.name));

probdata.transf_type = 3;
probdata.Ro_method   = 1;
probdata.flag_sens   = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DATA FIELDS IN 'ANALYSISOPT'  %%
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
analysisopt.i_max                = 500;      % Maximum number of iterations allowed in the search algorithm
analysisopt.e1                   = 0.001;    % Tolerance on how close design point is to limit-state surface
analysisopt.e2                   = 0.001;    % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code            = 0;        % 0: step size by Armijo rule, otherwise: given value is the step size
analysisopt.Recorded_u           = 1;        % 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
analysisopt.Recorded_x           = 1;        % 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations

% FORM, SORM analysis options
analysisopt.grad_flag            = 'ffd';    % 'ddm': direct differentiation, 'ffd': forward finite difference
analysisopt.ffdpara              = 1000;     % Parameter for computation of FFD estimates of gradients - Perturbation = stdv/analysisopt.ffdpara;
% Recommended values: 1000 for basic limit-state functions, 50 for FE-based limit-state functions
analysisopt.ffdpara_thetag       = 1000;     % Parameter for computation of FFD estimates of dbeta_dthetag
% perturbation = thetag/analysisopt.ffdpara_thetag if thetag ~= 0 or 1/analysisopt.ffdpara_thetag if thetag == 0;
% Recommended values: 1000 for basic limit-state functions, 100 for FE-based limit-state functions

% Simulation analysis (MC,IS,DS,SS) and distribution analysis options
analysisopt.num_sim              = 100000;   % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)

% Simulation analysis (MC, IS) and distribution analysis options
analysisopt.sim_point            = 'origin'; % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.0001;   % Target coefficient of variation for failure probability
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
%  DATA FIELDS IN 'GFUNDATA' (one structure per gfun)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator  = 'basic';
gfundata(1).type       = 'expression';   % Do no change this field!

% Expression of the limit-state function:
gfundata(1).expression = 'gfun_rigid_frame_system(M1, M2, M3, S1, S2, als)';

% Give explicit gradient expressions w.r.t. involved random variables (order in probdata.marg),
% if DDM is used (analysisopt.grad_flag='ddm')
% gfundata(1).dgdq       = { '  1 '
%                            ' -1 ' };

% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DATA FIELDS IN 'FEMODEL'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

femodel = [];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DATA FIELDS IN 'RANDOMFIELD'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomfield = [];

% =========================================================================
% RELI analysis
% =========================================================================
% analysisopt.form_solver             = 'interior-point';
% analysisopt.form_solver             = 'sqp';
analysisopt.form_solver             = 'cobyla';

% update before any analysis (must be run only once)
[probdata, gfundata, analysisopt]   = update_data(1, probdata, analysisopt, gfundata, []);

probdata.marg                       = distribution_parameter(probdata.marg);
probdata                            = det_rv_bounds(probdata);

n_ls                                = 6;
beta_1                              = nan(n_ls,1);
beta_3                              = nan(n_ls,1);
for ii = 1:n_ls    
    gfundata.cg                         = ii;
    formresults1                        = form(1, probdata, analysisopt, gfundata, [], []);
    
    formresults3                        = form_v3(1, probdata, analysisopt, gfundata, [], []);
    
    beta_1(ii) = formresults1.beta;
    beta_3(ii) = formresults3.beta;
    
end

for ii = 1:n_ls 
    disp(['limit state = ', num2str(ii)])
    disp(['beta_iRF = ', num2str(beta_1(ii)), ';     beta_', analysisopt.form_solver, ' = ', num2str(beta_3(ii))])
end