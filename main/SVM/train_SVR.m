function gfundata = train_SVR(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun

% Options settings

marg            = probdata.marg;
R               = probdata.correlation;
transf_type     = probdata.transf_type;
nrv             = size(marg,1);

NbCal           = analysisopt.NbCal;

flag_sens       = probdata.flag_sens;
sampling        = analysisopt.sampling;
block_size      = analysisopt.block_size;
rand_generator  = analysisopt.rand_generator;
% gridsel_option  = analysisopt.gridsel_option;
Sobolsetopt     = analysisopt.Sobolsetopt;
SVRbasis        = analysisopt.SVRbasis;
SVR_Nbasis      = analysisopt.SVR_Nbasis;
num_sim         = analysisopt.num_sim;
first_indices   = analysisopt.first_indices;
total_indices   = analysisopt.total_indices;
all_indices     = analysisopt.all_indices;

if isfield(gfundata(lsf),'ng')
   ng = gfundata(lsf).ng;
else
   ng = 1;
end

% Initialization

k  = 0;

S  = [];
ST = [];

Remp   = [];
indice = [];

% Sets the design of experiments for the learning basis

switch SVRbasis

   case 'Sobol_norm' % Points are deterministically generated (use of Sobol' sequences in hypercube [0; 1]) in the standard space (follow the normal distribution)

      if Sobolsetopt == 0 % Use of Statistical toolbox of Matlab with the sobolset function
         p  = sobolset(nrv);
         Ui = inv_norm_cdf(p(2:SVR_Nbasis+1,:)');
      end

      if Sobolsetopt == 1 % Use of Sobol' sequence
         a = zeros(nrv,SVR_Nbasis);
         for j = 1:SVR_Nbasis
            a(:,j) = sobolseq51(j,nrv)';
         end
         Ui = inv_norm_cdf(a);
      end
      
%       eval([ 'save -mat ' analysisopt.SVRbasis_filename ' Ui' ]);

   case 'CVT_unif' % Points are deterministicaly generated (use of Voronoï cells) in an hypersphere with a specified radius

      % Points are deterministicaly generated: use of Voronoi cells

      dim_num       = nrv;
      n             = SVR_Nbasis;
      batch         = 10000;
      init          = 3;
      sample        = 3;
      sample_num    = 10000;
      it_max        = 40;
      it_fixed      = 1;
%       init_string   = 'user';
%       sample_string = 'user';
      seed          = 123456789;
      r             = [];
      seed_init     = seed;
      [ r, seed, it_num, it_diff, energy ] = cvt( dim_num, n, batch, init, ...
         sample, sample_num, it_max, it_fixed, seed, r );

      % Choice of the radius of the hypersphere

      n_reg_radius = analysisopt.n_reg_radius;
      num = n_reg_radius;
      reg_radius = 0;
      block_u = inv_norm_cdf(twister(nrv,num));
      normuu = (dot(block_u,block_u)).^0.5;
      reg_radius = max(reg_radius,max(normuu));

      % Sample points in the standard u-space

      Ui = reg_radius * r;

   case 'file' % Ui points are recorded in a .mat file
      
      eval([ 'load ' analysisopt.SVRbasis_filename ]);

   otherwise
      
      error('train_SVR: analysisopt.SVRbasis should be one of ''Sobol_norm'' or ''Sobol_unif''')
      
end

% Model evaluation for the design of experiments

Xi = u_to_x(Ui,probdata);
[ allg, dummy ] = gfun(lsf,Xi,'no ',probdata,analysisopt,gfundata,femodel,randomfield);

for ig = 1:ng

   gridsel_option = 1;
   
   gmax = max(allg(ig,:));
   gmin = min(allg(ig,:));

   gfundata(lsf).SVR(ig).gmax = gmax;
   gfundata(lsf).SVR(ig).gmin = gmin;

   Gi   = ( allg(ig,:) - gmin ) / ( gmax - gmin );

   Di   = data(Ui',Gi'); % Learning basis


   % Cross validation for hyperparameters: C, epsilon and sigma (RBF is chosen for the kernel)
   % Range of parameters C, epsilon and sigma have to be defined by the user based on its knowledge of the problem

   for i = [ 1000000 ] % Choice for C [ 100 1000:1000:10000 10000:10000:100000 ]

      for j = -11:1:-1 % Choice for epsilon

         k = k + 1;

         rbf_tab = [ 1 2 5 7 8 10 15 20 30 50 ]; % Choice for sigma
         C  = i;
         epsilon = 2^j;

         SVR = svr_init(gridsel_option,rbf_tab,C,epsilon); % Select the best sigma (cross validation)

         CC(k) = i;
         eepsilon(k) = 2^j;

         [ r, SVR ] = train(SVR,Di); % Learning

         if gridsel_option == 1
            ind_SVR_best = SVR.best_index;
            indice(k)    = ind_SVR_best;
            Remp(1,k)    = SVR.best_score; % Minimal empirical error from cross validation
            SVR          = SVR.best;
         end

      end

   end

   % Select C, epsilon and sigma which minimizes the empirical error and from the cross validation

   [ remp_min, i_min ] = min(Remp);
   ind_SVR = indice(i_min);

   % Training with the best C, epsilon and sigma

   gridsel_option = 2;
   SVR = svr_init(gridsel_option,rbf_tab,CC(i_min),eepsilon(i_min),ind_SVR);
   [ r, SVR ] = train(SVR,Di);

   gfundata(lsf).SVR(ig).hyperparameters = SVR;
   gfundata(lsf).SVR(ig).error           = remp_min;

end

gfundata(lsf).evaluator           = 'svr';