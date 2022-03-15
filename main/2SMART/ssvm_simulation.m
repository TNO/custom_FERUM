function [ ssvmsimulationresults, probdata ] = ssvm_simulation(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun
nfun = 0;

% Number of random variables
nrv = size(probdata.marg,1); probdata.nrv = nrv;


% Cross-validation - Search of the best value for the sigma-parameter of the Gaussian kernel
gridsel_svc_flag = 1; % Set to 1 (do not modify)
rbf_tab0         = [ 1.5 3 5 7 10 15 20 30 50 100 200 400 ]; rbf_tab = rbf_tab0; ind_SVC0  = 0; % Range of values for sigma

% Sigma-parameter set to a specific value - Not used in 2SMART
% gridsel_svc_flag = 0;
% rbf_tab0 = [10]; rbf_tab = rbf_tab0;


analysisopt.n_max_for_clust = 10000;      % Set to 10000, do not modify
analysisopt.flag_cluster_supp_close = 0;  % Set to 0, do not modify

% Sampling strategy (percent values)
analysisopt.select_strat = [ ... 
% Stage :  Localization  |    Stabilization     |               Convergence
% Iter  :  1  n_supp1-1  |  n_supp1  n_supp2-1  |  n_supp2  n_supp3-1  n_supp3  n_supp_stop
         100        100          80         40           0          0        0            0; ...  % Cluster centers of work population points in the margin (maximum of 10000)
           0          0           0         10          10         10       10           10; ...  % Points in the margin which are the closest to the SVM classifier
           0          0          20         50          90         90       90           90 ];    % Cluster centers of points switching from one side of the classifier to another


analysisopt.flag_plotfigind      = 1; 
analysisopt.flag_plot23_all      = 1;
analysisopt.flag_plot23_lastonly = 1;

analysisopt.flag_pause           = 0; % 1: pause after each plot, 0: no pause
analysisopt.flag_mem             = 0; % Set to 0, do not modify

gca_param.auto_gca               = 0; % 1: automatic axis limits, 0: axis limits selected by user (variable AllLim in ssvm_subset_step.m)
gca_param.plot_thresholds        = 1; % 1: plot exact intermediate thresholds of the limit-state function for 2 or 3 dimensional random spaces
                                      %    Important: This is made by means of direct calls to the limit-state function (gfunbasic.m function) in plot_svm2d.m or plot_svm3d.m
                                      % 0: no plots
gca_param.plot_small             = 0;


if ~isfield(analysisopt,'flag_mem')
   analysisopt.flag_mem = 0;
end

if ~isfield(analysisopt,'flag_plot')
   analysisopt.flag_plot = 0;
end

if ~isfield(analysisopt,'flag_plot_gen')
   analysisopt.flag_plot_gen = 0;
end

if ~isfield(analysisopt,'flag_pause')
   analysisopt.flag_pause = 1;
end

if ~isfield(analysisopt,'echo_flag')
   analysisopt.echo_flag = 1;
end


rand_generator = analysisopt.rand_generator;

% Number of standard Gaussian samples for selection of Ru radius
% These are default values of the 2SMART algorithm - Please refer to Bourinet et al. paper for details
analysisopt.N_Ru_radius      = 200000;

% Numbers of samples for work populations
% These are default values of the 2SMART algorithm - Please refer to Bourinet et al. paper for details
analysisopt.nusvm1           = 10000;  % N1 parameter in Bourinet et al. paper
analysisopt.nusvm2           = 50000;  % N2 parameter in Bourinet et al. paper
analysisopt.nusvm3           = 200000; % N3 parameter in Bourinet et al. paper
analysisopt.nusvm4           = 200000; % N3 parameter in Bourinet et al. paper

% lambda-factors for lambdarj-mM algorithm
% These are default values of the 2SMART algorithm - Please refer to Bourinet et al. paper for details
analysisopt.ratio_amp_factor = [ 7 3.5 1 1 ];

% Number of iterations spent in the convergence stage
% These are default values of the 2SMART algorithm - Please refer to Bourinet et al. paper for details
analysisopt.N_supp_23        = [ 10 10 10 ];
analysisopt.N_supp_3stop     = [ 6 6 6 ];


allU   = [];
allG   = [];
allInd = [];
Nfun   = 0;

% Set central point for crude Monte Carlo
point = zeros(nrv,1); analysisopt.point = point;
   

% Modify correlation matrix and perform Cholesky decomposition
if ~isfield(probdata,'Lo')
    
   if probdata.transf_type == 3
    
      % Compute corrected correlation coefficients
      switch probdata.Ro_method
         case 0
            Ro = mod_corr( probdata.marg, probdata.correlation );
         case 1
            [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( probdata.marg, probdata.correlation , 0);
      end     
      probdata.Ro = Ro;
   
      % Cholesky decomposition
      % Lo = (chol(Ro))'; probdata.Lo = Lo;
      [Lo,ierr] = my_chol(Ro); if  ierr>0, return, end
      probdata.Lo = Lo;
   
      iLo = inv(Lo);
      probdata.iLo = iLo;

   end

end

% Establish covariance matrix, its Cholesky decomposition and its inverse
stdv1 = analysisopt.stdv_sim;
covariance = stdv1^2 * eye(nrv); probdata.covariance = covariance;
chol_covariance = stdv1 * eye(nrv);  probdata.chol_covariance = chol_covariance; % chol_covariance = chol(covariance);
inv_covariance = 1/stdv1 * eye(nrv); % inv_covariance = inv(covariance);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core part of the 2SMART algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First subset-like step ( i = 1 in Bourinet et al. paper )
Nb_step = 0;
ssvm_subset_step, eval(['save restart_from_step_' num2str(Nb_step) '.mat']); % pause

% Next subset-like steps ( i > 1 in Bourinet et al. paper )
% Loop until the limit-state function is equal to zero
while g0 > 0
   Nb_step = Nb_step + 1;
   ssvm_subset_step, eval(['save restart_from_step_' num2str(Nb_step) '.mat']); % pause
end


% Evaluation of the failure probability
pf     = 1;
pfbest = 1;
for istep = 0:Nb_step
   pf     = pf*svmdata(istep+1).svmpf(1,end);
   pfbest = pfbest*svmdata(istep+1).svmpf(1,svmdata(istep+1).bestIter);
end


% Store results
ssvmsimulationresults.pf      = pf;
ssvmsimulationresults.pfbest  = pfbest;
ssvmsimulationresults.beta    = -inv_norm_cdf(ssvmsimulationresults.pf);

ssvmsimulationresults.allU    = allU;
ssvmsimulationresults.allG    = allG;
ssvmsimulationresults.allInd  = allInd;
ssvmsimulationresults.nfun    = Nfun;
ssvmsimulationresults.Nb_step = Nb_step;

ssvmsimulationresults.svmdata = svmdata;
