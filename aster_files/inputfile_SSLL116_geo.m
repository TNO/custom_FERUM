% FORM results:
% 
% Number of iterations:
%      5
% 
% Reliability index beta:
%   2.6950e+000
% 
% Failure probability:
%   3.5189e-003
% 
% Number of calls to the limit-state function: 
%     35
% 
% formresults.dsptu
% 
% ans =
% 
%  -1.0060e+000
%   8.3370e-003
%  -1.6206e+000
%  -4.0119e-001
%  -1.8612e+000
% 
% formresults.dsptx
% 
% ans =
% 
%   1.7646e+011
%   3.0033e-001
%  -2.8103e+007
%   1.9198e-002
%   4.0694e-002
%   




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'PROBDATA'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.
probdata.name =  { 'h'
                   'E'
                   'NU'
                   'F_VERT'
                   'EP1'
                   'RAY1' };
              
% Marginal distributions for each random variable
% probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];
probdata.marg =    [ 0       0.1             0       0.1   nan   nan  nan  nan 0 ;
                     1  196200e6  196200e6*0.1  196200e6   nan   nan  nan  nan 0 ;
                     6       nan           nan       0.3  0.25  0.35  nan  nan 1 ;
                     1     -20e6     20e6*0.25     -20e6   nan   nan  nan  nan 0 ;
                     1      0.02      0.02*0.1      0.02   nan   nan  nan  nan 0 ;
                     1      0.05      0.05*0.1      0.05   nan   nan  nan  nan 0 ];
                   
% Correlation matrix
probdata.correlation = eye(6);

probdata.transf_type = 3;
probdata.Ro_method   = 1;
probdata.flag_sens   = 1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'ANALYSISOPT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'FEMODEL'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'GFUNDATA' (one structure per gfun)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'RANDOMFIELD'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomfield = [];