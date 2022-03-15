% Plotting options

flag_plotfigind      = analysisopt.flag_plotfigind; 
flag_plot23_all      = analysisopt.flag_plot23_all;
flag_plot23_lastonly = analysisopt.flag_plot23_lastonly;

if flag_plotfigind
   figure(100), clf
end

if flag_plot23_all | flag_plot23_lastonly
   if ~gca_param.auto_gca
      AllLim  = [-5 5];
      AllTick = [-5 -2.5 0 2.5 5];
      gca_param.XLim  = AllLim;
      gca_param.YLim  = AllLim;
      gca_param.ZLim  = AllLim;
      gca_param.XTick = AllTick;
      gca_param.YTick = AllTick;
      gca_param.ZTick = AllTick;
   end
end

nrv             = probdata.nrv;

flag_mem        = analysisopt.flag_mem;

rand_generator  = analysisopt.rand_generator;
width           = analysisopt.width;

nusvm1          = analysisopt.nusvm1;
nusvm2          = analysisopt.nusvm2;
nusvm3          = analysisopt.nusvm3;
nusvm4          = analysisopt.nusvm4;

n_max_for_clust = analysisopt.n_max_for_clust;

if isfield(analysisopt,'flag_var_rbf')
   flag_var_rbf = analysisopt.flag_var_rbf;
else
   flag_var_rbf = 0;
end

% Nu and Nn numbers of training samples - Please refer to Bourinet et al. paper for details
if Nb_step == 0
   num_sim_sdu  = analysisopt.num_sim_sdu(1);
   num_sim_norm = analysisopt.num_sim_norm(1);
else
   num_sim_sdu  = analysisopt.num_sim_sdu(2);
   num_sim_norm = analysisopt.num_sim_norm(2);
end


if Nb_step == 0

   % Find Ru radius of the hypersphere for further placements of uniform samples in it ( first subset-like step only - i = 1 )
   
   N_Ru_radius = analysisopt.N_Ru_radius;
   num_sim = N_Ru_radius;
   k = 0;
   Ru_radius = 0;
   block_size = analysisopt.buf_size;
   while k < num_sim
      block_size = min( block_size, num_sim-k );
      k = k + block_size;
      block_u = inv_norm_cdf(twister(nrv,block_size));
      normuu = (dot(block_u,block_u)).^0.5;
      Ru_radius = max(Ru_radius,max(normuu));
   end
   
   Ru_radius;

else

   SVC0 = svmdata(Nb_step).SVC;

end % if Nb_step == 0


% Selection of the Nu learning points uu - Please refer to Bourinet et al. paper for details
if Nb_step == 0
   
   % Uniform sampling in the hypersphere of given radius Ru

   num_sim = num_sim_sdu;

   k = 0;
   percent_done = 0;
   block_size = analysisopt.block_size;
   while k < num_sim

      block_size = min( block_size, num_sim-k );
      k = k + block_size;

      block_u = Ru_radius * ssvm_sdu(nrv,block_size);

      block_x = u_to_x(block_u,probdata);

      [ block_g, dummy ] = gfun(lsf,block_x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);

      allU = [ allU block_u ];
      allG = [ allG block_g ];

      if floor( k/num_sim * 20 ) > percent_done
         percent_done = floor( k/num_sim * 20 );
         fprintf(1,'Subset step #%d - %d%% complete\n',Nb_step,percent_done*5)
      end
      
   end

else
   
   germset = 1;

   % Load germs
   Subgerm0  = [];
   Subgermg0 = [];
   for ik = 1:Ikmax(germset)
      fprintf(1,[ 'Germ from file step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat\n' ]);
      eval([ 'load step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat' ]);
      Subgerm0  = [ Subgerm0 subgerm0 ];
      Subgermg0 = [ Subgermg0 subgermg0 ];
   end

   % Generation of num_sim = nusvm1 points using the lambda rj-mM algorithm based on the current SVM classifier (please refer to Bourinet et al. Paper)
   num_sim = nusvm1;
   ratio_amp_factor = analysisopt.ratio_amp_factor(germset); % lambda factor = 7
   genpopsubset

   % Determination of the Nu learning points uu from these num_sim = nusvm1 points by clustering (kmeans algorithmm)
   
   k = 0;
   ik = 0;
   subset = [];
   block_size = analysisopt.buf_size;
   while k < size_pop_subset
      block_size = min( block_size, size_pop_subset-k );
      k = k + block_size;
      ik = ik + 1;
      fprintf(1,[ 'load step' num2str(Nb_step) '_pop_' num2str(ik) '.mat\n' ]);
      eval([ 'load step' num2str(Nb_step) '_pop_' num2str(ik) '.mat' ]);
      subset = [ subset block_u ];      
   end
   
   num_sim = num_sim_sdu;
   
   fprintf(1,'\nStart clustering stage\n');
   tic;
   [ centers, mincenter, mindist, q2, quality ] = ucsd_kmeans(subset',num_sim+1);
   t = toc; fprintf(1,'Time required for clustering %d SVM subset points into %d computation points (sec): %f\n\n', size(subset,2), num_sim, t);
   subsubset = centers(2:end,:)';

   % Call the limit-state function for the Nu learning points uu

   num_sim = num_sim_sdu;
   
   k = 0;
   percent_done = 0;
   block_size = analysisopt.block_size;
   while k < num_sim
      
      block_size = min( block_size, num_sim-k );
      k = k + block_size;
      
      block_u = subsubset(:,(k-block_size+1):k);
      
      block_x = u_to_x(block_u,probdata);
      
      [ block_g, dummy ] = gfun(lsf,block_x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
      
      allU = [ allU block_u ];
      allG = [ allG block_g ];

      if floor( k/num_sim * 20 ) > percent_done
         percent_done = floor( k/num_sim * 20 );
         fprintf(1,'Subset step #%d - %d%% complete\n',Nb_step,percent_done*5)
      end

   end

end % if Nb_step == 0

allInd = [ allInd; Nfun+1 Nfun+num_sim ];
Nfun = Nfun + num_sim;


% Selection of the Nn learning points un - Please refer to Bourinet et al. paper for details
if Nb_step > 0

   germset = 1;

   % Load germs
   Subgerm0 = [];
   Subgermg0 = [];
   for ik = 1:Ikmax(germset)
      fprintf(1,[ 'Germ from file step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat\n' ]);
      eval([ 'load step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat' ]);
      Subgerm0 = [ Subgerm0 subgerm0 ];
      Subgermg0 = [ Subgermg0 subgermg0 ];
   end

   % Generation of num_sim = nusvm1 points using the lambda rj-mM algorithm based on the current SVM classifier (please refer to Bourinet et al. Paper)
   num_sim = nusvm1;
   ratio_amp_factor = 1; % lambda factor = 1 - Original Modified Metropolis Hasting algorithm (termed rj-mM algorithm in Bourinet et al. paper)
   genpopsubset

   % Determination of the Nn learning points un from these num_sim = nusvm1 points (random selection)
   
   k =  0;
   ik = 0;
   subset = [];
   block_size = analysisopt.buf_size;
   while k < size_pop_subset
      block_size = min( block_size, size_pop_subset-k );
      k = k + block_size;
      ik = ik + 1;
      fprintf(1,[ 'load step' num2str(Nb_step) '_pop_' num2str(ik) '.mat\n' ]);
      eval([ 'load step' num2str(Nb_step) '_pop_' num2str(ik) '.mat' ]);
      subset = [ subset block_u ];
   end

   Inorm = randperm(size_pop_subset);
   subsubset = subset(:,Inorm(1:num_sim_norm));

end % if Nb_step > 0


% Call the limit-state function for the Nn learning points un

num_sim = num_sim_norm;

k = 0;
percent_done = 0;
block_size = analysisopt.block_size;
while k < num_sim

   block_size = min( block_size, num_sim-k );
   k = k + block_size;

   if Nb_step == 0
      block_u = inv_norm_cdf(twister(nrv,block_size));
   else
      block_u = subsubset(:,(k-block_size+1):k);
   end

   block_x = u_to_x(block_u,probdata);

   [ block_g, dummy ] = gfun(lsf,block_x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);

   allU = [ allU block_u ];
   allG = [ allG block_g ];

   if floor( k/num_sim * 20 ) > percent_done
      percent_done = floor( k/num_sim * 20 );
      fprintf(1,'Subset step #%d - %d%% complete\n',Nb_step,percent_done*5)
   end

end

allInd = [ allInd; Nfun+1 Nfun+num_sim ];
Nfun = Nfun + num_sim;


pf_target = analysisopt.pf_target;

% Find the limit-state threshold value yi as the alpha-quantile on the Nn learning points un
g0  = ferum_quantile(allG((Nfun-num_sim_norm+1):Nfun),pf_target);
pf0 = pf_target;
if g0 < 0
   g0  = 0;
   pf0 = length( find( allG((Nfun-num_sim_norm+1):Nfun) < g0 ) ) / num_sim_norm;
end
fprintf(1,'\n\n>>> g0 = %f  -  pf0 = %f\n', g0, pf0);

if g0 > 0
   g0_try  = ferum_quantile(allG((Nfun-num_sim_norm+1):Nfun),pf_target*3/4);
   pf0_try = pf_target*3/4;
   if g0_try < 0
      g0_try  = 0;
      pf0_try = length( find( allG((Nfun-num_sim_norm+1):Nfun) < g0_try ) ) / num_sim_norm;
   end
   if g0_try == 0
      g0  = g0_try;
      pf0 = pf0_try;
      fprintf(1,'>>> g0 = %f  -  pf0 = %f\n', g0, pf0);
   end
end
fprintf(1,'\n\n');
         

if Nb_step == 0

   I_train = 1:length(allG);

else

   gridsel_option_svc = 0;
   num_sim = length(allG);
   k = 0;
   I_train = [];
   block_size = analysisopt.buf_size;
   svm_buf_size = analysisopt.svm_buf_size;
   while k < num_sim
      block_size = min( block_size, num_sim-k );
      k = k + block_size;
      block_u = allU(:,(k-block_size+1):k);
      block_g = eval_svm(block_u,svm_buf_size,SVC0,gridsel_option_svc);
      I_train = [ I_train find( block_g < 0 ) ];
   end

end

% First SVM classifier based on the Nu and Nn learning points
if Nb_step == 0
   if gridsel_svc_flag
      gridsel_option_svc = 1;
   else
      gridsel_option_svc = 0;
   end
   rbf_tab = rbf_tab0;
else
   if flag_var_rbf == 1
      if gridsel_svc_flag
         gridsel_option_svc = 0;
         rbf_tab = rbf_tab0(ind_SVC0+ind_SVC_best);
      else
         gridsel_option_svc = 0;
      end
   elseif flag_var_rbf == 0
      if gridsel_svc_flag
         gridsel_option_svc = 0;
         rbf_tab = rbf_tab0(ind_SVC0+ind_SVC_best);
      else
         gridsel_option_svc = 0;
      end
   end
end
SVC = svm_init(gridsel_option_svc,rbf_tab);
Y = sign( allG(1,I_train) - g0 );
D = data(allU(:,I_train)',Y')
[ r, SVC ] = train(SVC,D);
if gridsel_option_svc == 1
   ind_SVC_best = SVC.best_index;
   gridsel_option_svc = 2;
   SVC = svm_init(gridsel_option_svc,rbf_tab,ind_SVC_best);
   [ r, SVC ] = train(SVC,D);
end


% Initialize 
select_strat            = analysisopt.select_strat;
flag_cluster_supp_close = analysisopt.flag_cluster_supp_close;

if Nb_step == 0
   n_supp_23    = analysisopt.N_supp_23(1);
   n_supp_3stop = analysisopt.N_supp_3stop(1);
else
   if g0 > 0
      n_supp_23    = analysisopt.N_supp_23(2);
      n_supp_3stop = analysisopt.N_supp_3stop(2);
   else
      n_supp_23    = analysisopt.N_supp_23(3);
      n_supp_3stop = analysisopt.N_supp_3stop(3);
   end
end
n_supp2     = round(12*(nrv/2)^0.2);
n_supp1     = round(6*(nrv/2)^0.2);
n_supp3     = n_supp2 + n_supp_23;
n_supp_stop = n_supp3 + n_supp_3stop;

fprintf(1,[ '\nn_supp1 = ' num2str(n_supp1) ' - n_supp2 = ' num2str(n_supp2) ' - n_supp3 = ' num2str(n_supp3) ' - n_supp_stop = ' num2str(n_supp_stop) '\n' ]);

if nrv==2 & flag_plot23_all

   figure(20+Nb_step)
   clf
   plot_svm2d(SVC,gca_param);

   if gca_param.auto_gca
      set(gca,'FontSize',14);
      axis equal
   else
      set(gca,'FontSize',14,'XLim',gca_param.XLim,'YLim',gca_param.YLim,'XTick',gca_param.XTick,'YTick',gca_param.YTick);
      axis square
   end
   grid on
   box on

   xlabel(probdata.name(1),'FontSize',14);
   ylabel(probdata.name(2),'FontSize',14);

   if analysisopt.flag_pause
      disp('Press a key'), pause
   else
      pause(1)
   end

elseif nrv==3 & flag_plot23_all

   figure(20+Nb_step)
   clf
   plot_svm3d(SVC,gca_param);

   if gca_param.auto_gca
      set(gca,'FontSize',14);
      axis equal
   else
      set(gca,'FontSize',14,'XLim',gca_param.XLim,'YLim',gca_param.YLim,'ZLim',gca_param.ZLim,'XTick',gca_param.XTick,'YTick',gca_param.YTick,'ZTick',gca_param.ZTick);
      axis square
   end
   grid on
   box on

   xlabel(probdata.name(1),'FontSize',14);
   ylabel(probdata.name(2),'FontSize',14);
   zlabel(probdata.name(3),'FontSize',14);

   if analysisopt.flag_pause
      disp('Press a key'), pause
   else
      pause(1)
   end

end

% Main loop - k-iterations in Bourinet et al. paper
kk               = 0;
svmpf            = nan*zeros(3,n_supp_stop+1);
svmmargin        = nan*zeros(2,n_supp_stop);
svmswitch        = nan*zeros(2,n_supp_stop);
svmmisclassified = nan*zeros(4,n_supp_stop);
while kk < n_supp_stop

   kk = kk + 1;
   fprintf(1,'\nIteration : %d\n', kk);

   if kk <= n_supp1
      nusvm = nusvm1;
      germset = 1;
      n_pt_supp        = round( ( 3.0 + 1.0*(kk-1)/(n_supp1-1) ) * nrv^0.5 );
      n_pt_supp_switch = ceil( ( select_strat(3,1) + (kk-1)/(n_supp1-1) * (select_strat(3,2)-select_strat(3,1)) ) /100 * n_pt_supp );
      if mod(kk-0,2) == 0
         n_pt_supp_clust  = ceil( ( select_strat(1,1) + (kk-1)/(n_supp1-1) * (select_strat(1,2)-select_strat(1,1)) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,1) + (kk-1)/(n_supp1-1) * (select_strat(2,2)-select_strat(2,1)) ) /100 * n_pt_supp );
      else
         n_pt_supp_switch = 0;
         n_pt_supp_clust  = ceil( ( select_strat(1,1)/(select_strat(1,1)+select_strat(2,1))*100 + (kk-1)/(n_supp1-1) * (select_strat(1,2)/(select_strat(1,2)+select_strat(2,2))*100-select_strat(1,1)/(select_strat(1,1)+select_strat(2,1))*100) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,1)/(select_strat(1,1)+select_strat(2,1))*100 + (kk-1)/(n_supp1-1) * (select_strat(2,2)/(select_strat(1,2)+select_strat(2,2))*100-select_strat(2,1)/(select_strat(1,1)+select_strat(2,1))*100) ) /100 * n_pt_supp );
      end
   elseif kk <= n_supp2
      nusvm = nusvm2;
      germset = 2;
      n_pt_supp        = round( ( 4.0 + 1.0*(kk-n_supp1)/(n_supp2-n_supp1) ) * nrv^0.5 );
      n_pt_supp_switch = ceil( ( select_strat(3,3) + (kk-n_supp1)/(n_supp2-n_supp1) * (select_strat(3,4)-select_strat(3,3)) ) /100 * n_pt_supp );
      if mod(kk-n_supp1,2) == 0
         n_pt_supp_clust  = ceil( ( select_strat(1,3) + (kk-n_supp1)/(n_supp2-n_supp1) * (select_strat(1,4)-select_strat(1,3)) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,3) + (kk-n_supp1)/(n_supp2-n_supp1) * (select_strat(2,4)-select_strat(2,3)) ) /100 * n_pt_supp );
      else
         n_pt_supp_switch = 0;
         n_pt_supp_clust  = ceil( ( select_strat(1,3)/(select_strat(1,3)+select_strat(2,3))*100 + (kk-n_supp1)/(n_supp2-n_supp1) * (select_strat(1,4)/(select_strat(1,4)+select_strat(2,4))*100-select_strat(1,3)/(select_strat(1,3)+select_strat(2,3))*100) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,3)/(select_strat(1,3)+select_strat(2,3))*100 + (kk-n_supp1)/(n_supp2-n_supp1) * (select_strat(2,4)/(select_strat(1,4)+select_strat(2,4))*100-select_strat(2,3)/(select_strat(1,3)+select_strat(2,3))*100) ) /100 * n_pt_supp );
      end
   elseif kk <= n_supp3
      nusvm = nusvm3;
      germset = 3;
      n_pt_supp        = round( 5.0 * nrv^0.5 );
      n_pt_supp_switch = ceil( ( select_strat(3,5) + (kk-n_supp2)/(n_supp3-n_supp2) * (select_strat(3,6)-select_strat(3,5)) ) /100 * n_pt_supp );
      if mod(kk-n_supp2,2) == 0
         n_pt_supp_clust  = ceil( ( select_strat(1,5) + (kk-n_supp2)/(n_supp3-n_supp2) * (select_strat(1,6)-select_strat(1,5)) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,5) + (kk-n_supp2)/(n_supp3-n_supp2) * (select_strat(2,6)-select_strat(2,5)) ) /100 * n_pt_supp );
      else
         n_pt_supp_switch = 0;
         n_pt_supp_clust  = ceil( ( select_strat(1,5)/(select_strat(1,5)+select_strat(2,5))*100 + (kk-n_supp2)/(n_supp3-n_supp2) * (select_strat(1,6)/(select_strat(1,6)+select_strat(2,6))*100-select_strat(1,5)/(select_strat(1,5)+select_strat(2,5))*100) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,5)/(select_strat(1,5)+select_strat(2,5))*100 + (kk-n_supp2)/(n_supp3-n_supp2) * (select_strat(2,6)/(select_strat(1,6)+select_strat(2,6))*100-select_strat(2,5)/(select_strat(1,5)+select_strat(2,5))*100) ) /100 * n_pt_supp );
      end
   else
      nusvm = nusvm4;
      if nusvm4 == nusvm3
         germset = 3;
      else
         germset = 4;
      end
      n_pt_supp        = round( 5.0 * nrv^0.5 );
      n_pt_supp_switch = ceil( ( select_strat(3,7) + (kk-n_supp3)/(n_supp_stop-n_supp3) * (select_strat(3,8)-select_strat(3,7)) ) /100 * n_pt_supp );
      if mod(kk-n_supp3,2) == 0
         n_pt_supp_clust  = ceil( ( select_strat(1,7) + (kk-n_supp3)/(n_supp_stop-n_supp3) * (select_strat(1,8)-select_strat(1,7)) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,7) + (kk-n_supp3)/(n_supp_stop-n_supp3) * (select_strat(2,8)-select_strat(2,7)) ) /100 * n_pt_supp );
      else
         n_pt_supp_switch = 0;
         n_pt_supp_clust  = ceil( ( select_strat(1,7)/(select_strat(1,7)+select_strat(2,7))*100 + (kk-n_supp3)/(n_supp_stop-n_supp3) * (select_strat(1,8)/(select_strat(1,8)+select_strat(2,8))*100-select_strat(1,7)/(select_strat(1,7)+select_strat(2,7))*100) ) /100 * n_pt_supp );
         n_pt_supp_close  = ceil( ( select_strat(2,7)/(select_strat(1,7)+select_strat(2,7))*100 + (kk-n_supp3)/(n_supp_stop-n_supp3) * (select_strat(2,8)/(select_strat(1,8)+select_strat(2,8))*100-select_strat(2,7)/(select_strat(1,7)+select_strat(2,7))*100) ) /100 * n_pt_supp );
      end
   end

   if Nb_step > 0

      % Load germs
      Subgerm0  = [];
      Subgermg0 = [];
      for ik = 1:Ikmax(germset)
         fprintf(1,[ 'Germ from file step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat\n' ]);
         eval([ 'load step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat' ]);
         Subgerm0  = [ Subgerm0 subgerm0 ];
         Subgermg0 = [ Subgermg0 subgermg0 ];
      end

     %  Generation of num_sim = nusvm points(N1, N2 or N3) using the lambda rj-mM algorithm based on the current SVM classifier (please refer to Bourinet et al. Paper)
      num_sim = nusvm;
      ratio_amp_factor = analysisopt.ratio_amp_factor(germset); % lambda factor = 7, 3.5 or 1 respectively for N1, N2 or N3
      genpopsubset

   end

   % Evaluation of num_sim points based on the SVM classifier
   if Nb_step == 0
      num_sim = nusvm;
   else
      num_sim = size_pop_subset;
   end

   if flag_mem
      Usvm           = zeros(nrv,num_sim);
      Usvm_in_margin = [];
      Usvm_switch    = [];
   else
      Usvm           = [];
      Usvm_in_margin = [];
      Usvm_switch    = [];
   end
   Gsvm              = zeros(1,num_sim);
   Gsvm_in_margin    = [];
   Gsvm_switch       = [];

   k                 = 0;
   ik                = 0;
   n_in_margin       = 0;
   n_switch          = 0;
   nusvm_reg_radius  = 0;
   block_size        = analysisopt.buf_size;
   svm_buf_size      = analysisopt.svm_buf_size;
   while k < num_sim

      block_size = min( block_size, num_sim-k );
      k = k + block_size;
      ik = ik + 1;

      if n_pt_supp_switch > 0
         eval(['load svm' num2str(ik) '.mat']);
         I_neg_previous = I_neg;
      else
         if Nb_step == 0
            if kk <= n_supp2
               block_u = Ru_radius * ssvm_sdu(nrv,block_size);
            else
               block_u = inv_norm_cdf(twister(nrv,block_size));
            end
         else
            eval([ 'load step' num2str(Nb_step) '_pop_' num2str(ik) '.mat' ]);
         end
      end
      block_g = eval_svm(block_u,svm_buf_size,SVC,gridsel_option_svc);

      if flag_mem
         Usvm(:,(k-block_size+1):k) = block_u;
      end
      Gsvm((k-block_size+1):k) = block_g;

      normuu = (dot(block_u,block_u)).^0.5;
      nusvm_reg_radius = max(nusvm_reg_radius , max(normuu));

      % Points in the margin
      I_in_margin = find( 1 > abs(block_g) );
      n_in_margin = n_in_margin + length(I_in_margin);
      if flag_mem
         Usvm_in_margin = [ Usvm_in_margin block_u(:,I_in_margin) ];
      end
      Gsvm_in_margin = [ Gsvm_in_margin block_g(I_in_margin) ];

      I_neg = find( block_g < 0 );

      % Switching points
      if n_pt_supp_switch > 0
         test_neg_previous = zeros(1,block_size); test_neg_previous(I_neg_previous) = 1;
         test_neg = zeros(1,block_size); test_neg(I_neg) = 1;
         I_switch = find( ( test_neg_previous + test_neg ) == 1 );
         n_switch = n_switch + length(I_switch);
         if flag_mem
            Usvm_switch = [ Usvm_switch block_u(:,I_switch) ];
         end
         Gsvm_switch = [ Gsvm_switch block_g(I_switch) ];
      end

      if flag_mem
         eval(['save svm' num2str(ik) '.mat block_u I_neg']);
      else
         if n_pt_supp_switch > 0
            eval(['save svm' num2str(ik) '.mat block_u block_g I_in_margin I_neg I_switch']);
         else
            eval(['save svm' num2str(ik) '.mat block_u block_g I_in_margin I_neg']);
         end
      end

   end
   ikmax = ik;

   n_def_U      = length(find(Gsvm<=0));  % Number of failure points
   n_inf_minus1 = length(find(Gsvm<=-1)); % Number of points < -1
   n_inf_plus1  = length(find(Gsvm<=1));  % Number of points < +1
   n_U          = length(Gsvm);           % Total number of points (this number is different from nusvm)
   fprintf(1,'\nNumber of failure points             : %d\n', n_def_U);
   fprintf(1,'Total number of learning points      : %d\n', n_U)
   fprintf(1,'Maximal radius of learning points    : %e\n', nusvm_reg_radius);
   if kk > n_supp2
      fprintf(1,'Pf estimate                          : %e\n', n_def_U/n_U);
      fprintf(1,'Pf < [ -1 , 0 , +1 ]                 : [ %e , %e , %e ]\n', n_inf_minus1/n_U,n_def_U/n_U,n_inf_plus1/n_U);
      svmpf(1,kk) = n_def_U/n_U;
      svmpf(2,kk) = n_inf_minus1/n_U;
      svmpf(3,kk) = n_inf_plus1/n_U;
   end

   fprintf(1,'Number of points in the margin       : %d\n',n_in_margin);
   fprintf(1,'Ratio of points in the margin / total: %f\n',n_in_margin/n_U);
   svmmargin(1,kk) = n_in_margin;
   svmmargin(2,kk) = n_in_margin/n_U;

   if n_pt_supp_switch > 0
      fprintf(1,'Number of switching points           : %d\n',n_switch);
      fprintf(1,'Ratio of switching points / total    : %f\n',n_switch/n_U);
      svmswitch(1,kk) = n_switch;
      svmswitch(2,kk) = n_switch/n_U;
   end

   if kk > n_supp2
      if flag_plotfigind
         figure(100), subplot(2,2,1)
         hold on
         plot(1:kk,svmpf(1,1:kk),'ko-')
         plot(1:kk,svmpf(2,1:kk),'bo-')
         plot(1:kk,svmpf(3,1:kk),'ro-')
         xlabel('Iteration k');
         ylabel('pf');
         figure(100), subplot(2,2,3)
         hold on
         plot(1:kk,[ nan*ones(1,n_supp2) svmmargin(2,(n_supp2+1):kk) ],'ko-')
         xlabel('Iteration k');
         ylabel('Ratio margin / total');
         figure(100), subplot(2,2,4)
         hold on
         plot(1:kk,[ nan*ones(1,n_supp2) svmswitch(2,(n_supp2+1):kk) ],'ko-')
         xlabel('Iteration k');
         ylabel('Ratio switching / total');
      end
      if kk > (n_supp2+1)
         if svmmargin(2,kk) < min(svmmargin(2,(n_supp2+1):(kk-1)))
            bestSVC = SVC;
            bestIter = kk;
         end
      end
   end

   U = [];

   if flag_mem
      % Suppress redundant points in the margin, if any
      Isame = find( all( abs( Usvm_in_margin(:,2:n_in_margin) - Usvm_in_margin(:,1:n_in_margin-1) ) < eps ) );
      if ~isempty(Isame)
         Usvm_in_margin(:,Isame) = [];
         Gsvm_in_margin(Isame) = [];
      end
      n_in_margin_tmp = size(Usvm_in_margin,2);
   end
   
   if n_pt_supp_clust > 0
      % Clustering of points in the margin -> n_pt_supp_clust points
      fprintf(1,'\nStart clustering stage\n'); tic;
      if flag_mem
      else
         Usvm_in_margin = [];
         for ik = 1:ikmax
            eval(['load svm' num2str(ik) '.mat']);
            Usvm_in_margin = [ Usvm_in_margin block_u(:,I_in_margin) ];
            n_in_margin_tmp = size(Usvm_in_margin,2);
            Isame = find( all( abs( Usvm_in_margin(:,2:n_in_margin_tmp) - Usvm_in_margin(:,1:n_in_margin_tmp-1) ) < eps ) );
            if ~isempty(Isame)
               Usvm_in_margin(:,Isame) = [];
               n_in_margin_tmp = size(Usvm_in_margin,2);
            end
            if n_in_margin_tmp >= n_max_for_clust, break, end
         end
      end
      if n_in_margin_tmp > n_max_for_clust
         [ centers, mincenter, mindist, q2, quality ] = ucsd_kmeans(Usvm_in_margin(:,1:n_max_for_clust)',n_pt_supp_clust+1);
         U_clust = centers(2:end,:)';
         for jc = 1:n_pt_supp_clust
            norm_clust = dot( (U_clust(:,jc)*ones(1,n_max_for_clust)-Usvm_in_margin(:,1:n_max_for_clust)) , ...
               (U_clust(:,jc)*ones(1,n_max_for_clust)-Usvm_in_margin(:,1:n_max_for_clust)) );
            [dummy,j_clust] = min(norm_clust);
            U_clust(:,jc) = Usvm_in_margin(:,j_clust);
         end
      elseif n_in_margin_tmp > n_pt_supp_clust
         [ centers, mincenter, mindist, q2, quality ] = ucsd_kmeans(Usvm_in_margin(:,1:n_in_margin_tmp)',n_pt_supp_clust+1);
         U_clust = centers(2:end,:)';
         for jc = 1:n_pt_supp_clust
            norm_clust = dot( (U_clust(:,jc)*ones(1,n_in_margin_tmp)-Usvm_in_margin(:,1:n_in_margin_tmp)) , ...
               (U_clust(:,jc)*ones(1,n_in_margin_tmp)-Usvm_in_margin(:,1:n_in_margin_tmp)) );
            [dummy,j_clust] = min(norm_clust);
            U_clust(:,jc) = Usvm_in_margin(:,j_clust);
         end
      else
         U_clust = Usvm_in_margin;
         n_pt_supp_clust = n_in_margin_tmp;
      end
      U = [U U_clust];
      t = toc; fprintf(1,'Time required for clustering %d SVM subset points into %d computation points (sec): %f\n\n', min([n_max_for_clust n_in_margin_tmp]), n_pt_supp_clust, t);
   end % if n_pt_supp_clust > 0
   
   if n_pt_supp_close > 0
      % Points in the margin the closest to the classifier -> n_pt_supp_close points
      if flag_cluster_supp_close
         times_supp_close = 5;
      else
         times_supp_close = 1;
      end
      if flag_mem
         [ dummy, index ] = sort(abs(Gsvm_in_margin));
         Usvm_in_margin_sorted = Usvm_in_margin(:,index);
      else
         Usvm_in_margin_sorted = [];
         Gsvm_in_margin_sorted = [];
         for ik = 1:ikmax
            eval(['load svm' num2str(ik) '.mat']);
            Usvm_in_margin_sorted = [ Usvm_in_margin_sorted block_u(:,I_in_margin) ];
            Gsvm_in_margin_sorted = [ Gsvm_in_margin_sorted block_g(I_in_margin) ];
            n_in_margin_tmp = size(Usvm_in_margin_sorted,2);
            Isame = find( all( abs( Usvm_in_margin_sorted(:,2:n_in_margin_tmp) - Usvm_in_margin_sorted(:,1:n_in_margin_tmp-1) ) < eps ) );
            if ~isempty(Isame)
               Usvm_in_margin_sorted(:,Isame) = [];
               Gsvm_in_margin_sorted(Isame)   = [];
               n_in_margin_tmp = size(Usvm_in_margin_sorted,2);
            end
            [ dummy, index ] = sort(abs(Gsvm_in_margin_sorted));
            if n_in_margin_tmp > times_supp_close * n_pt_supp_close
               Usvm_in_margin_sorted = Usvm_in_margin_sorted(:,index(1:times_supp_close*n_pt_supp_close));
               Gsvm_in_margin_sorted = Gsvm_in_margin_sorted(index(1:times_supp_close*n_pt_supp_close));
            else
               Usvm_in_margin_sorted = Usvm_in_margin_sorted(:,index(1:n_in_margin_tmp));
               Gsvm_in_margin_sorted = Gsvm_in_margin_sorted(index(1:n_in_margin_tmp));
            end
            n_in_margin_tmp = size(Usvm_in_margin_sorted,2);
         end
      end
      if n_in_margin_tmp > n_pt_supp_close
         if flag_cluster_supp_close
            [ centers, mincenter, mindist, q2, quality ] = ucsd_kmeans(Usvm_in_margin_sorted(:,1:n_in_margin_tmp)',n_pt_supp_close+1);
            U_close = centers(2:end,:)';
            U_close0 = U_close;
            for jc = 1:n_pt_supp_close
               norm_clust = dot( (U_close(:,jc)*ones(1,n_in_margin_tmp)-Usvm_in_margin_sorted(:,1:n_in_margin_tmp)) , ...
                  (U_close(:,jc)*ones(1,n_in_margin_tmp)-Usvm_in_margin_sorted(:,1:n_in_margin_tmp)) );
               [dummy,j_close] = min(norm_clust);
               U_close(:,jc) = Usvm_in_margin_sorted(:,j_close);
            end
            Isame = find( all( abs( U_close(:,2:n_pt_supp_close) - U_close(:,1:n_pt_supp_close-1) ) < eps ) );
            if ~isempty(Isame)
               U_close(:,Isame) = [];
               n_pt_supp_close
               n_pt_supp_close = size(U_close,2);
               n_pt_supp_close
               disp('n_pt_supp_close : same points\n'); pause
            end
         else
            U_close = Usvm_in_margin_sorted(:,1:n_pt_supp_close);
         end
      else
         U_close = Usvm_in_margin_sorted(:,1:n_in_margin_tmp);
         n_pt_supp_close = n_in_margin_tmp;
      end
      U = [U U_close];
   end % if n_pt_supp_close > 0
   
   if n_pt_supp_switch > 0
      % Clustering of switching points -> n_pt_supp_switch points
      fprintf(1,'Start clustering stage\n'); tic;
      if flag_mem
         n_switch_tmp = n_switch;
      else
         Usvm_switch = [];
         for ik = 1:ikmax
            eval(['load svm' num2str(ik) '.mat']);
            Usvm_switch = [ Usvm_switch block_u(:,I_switch) ];
            n_switch_tmp = size(Usvm_switch,2);
            Isame = find( all( abs( Usvm_switch(:,2:n_switch_tmp) - Usvm_switch(:,1:n_switch_tmp-1) ) < eps ) );
            if ~isempty(Isame)
               Usvm_switch(:,Isame) = [];
               n_switch_tmp = size(Usvm_switch,2);
            end
            if n_switch_tmp >= n_max_for_clust, break, end
         end
      end
      if n_switch_tmp > n_max_for_clust
         [ centers, mincenter, mindist, q2, quality ] = ucsd_kmeans(Usvm_switch(:,1:n_max_for_clust)',n_pt_supp_switch+1);
         U_switch = centers(2:end,:)';
         for js = 1:n_pt_supp_switch
            norm_switch = dot( (U_switch(:,js)*ones(1,n_max_for_clust)-Usvm_switch(:,1:n_max_for_clust)) , ...
               (U_switch(:,js)*ones(1,n_max_for_clust)-Usvm_switch(:,1:n_max_for_clust)) );
            [dummy,j_switch] = min(norm_switch);
            U_switch(:,js) = Usvm_switch(:,j_switch);
         end
      elseif n_switch_tmp > n_pt_supp_switch
         [ centers, mincenter, mindist, q2, quality ] = ucsd_kmeans(Usvm_switch(:,1:n_switch_tmp)',n_pt_supp_switch+1);
         U_switch = centers(2:end,:)';
         for js = 1:n_pt_supp_switch
            norm_switch = dot( (U_switch(:,js)*ones(1,n_switch_tmp)-Usvm_switch(:,1:n_switch_tmp)) , ...
               (U_switch(:,js)*ones(1,n_switch_tmp)-Usvm_switch(:,1:n_switch_tmp)) );
            [dummy,j_switch] = min(norm_switch);
            U_switch(:,js) = Usvm_switch(:,j_switch);
         end
      else
         U_switch = Usvm_switch;
         n_pt_supp_switch = n_switch_tmp;
      end
      U = [U U_switch];
      t = toc; fprintf(1,'Time required for clustering %d SVM subset points into %d computation points (sec): %f\n\n', min([n_max_for_clust n_switch_tmp]), n_pt_supp_switch, t);
   end % if n_pt_supp_switch > 0

   if flag_plotfigind
      figure(100), subplot(2,2,2)
      hold on
      plot(kk,n_pt_supp,'k+')
      plot(kk,n_pt_supp_clust+n_pt_supp_close+n_pt_supp_switch,'k*')
      plot(kk,n_pt_supp_clust,'ro')
      plot(kk,n_pt_supp_close,'go')
      plot(kk,n_pt_supp_switch,'bo')
      xlabel('Iteration k');
      ylabel('Nb of learning points');
end

   if nrv==2 & flag_plot23_all

      figure(20+Nb_step)
      clf
      if gca_param.plot_thresholds
         for istep = 0:(Nb_step-1)
            G0(istep+1) = svmdata(istep+1).g0;
         end
         G0(Nb_step+1) = g0;
         if n_pt_supp_clust > 0
            plot_svm2d(SVC,gca_param, ...
                       Usvm,n_pt_supp_clust,Usvm_in_margin,[],n_pt_supp_switch,Usvm_switch, ...
                       probdata,analysisopt,gfundata,Nb_step,G0);
         else
            plot_svm2d(SVC,gca_param, ...
                       Usvm,n_pt_supp_clust,[],Usvm_in_margin_sorted,n_pt_supp_switch,Usvm_switch, ...
                       probdata,analysisopt,gfundata,Nb_step,G0);
         end
      else
         plot_svm2d(SVC,gca_param, ...
                    Usvm,n_pt_supp_clust,Usvm_in_margin,Usvm_in_margin_sorted,n_pt_supp_switch,Usvm_switch);
      end

      for iii = 1:size(allInd,1)
         if ( allInd(iii,2)-allInd(iii,1)+1 ) == num_sim_sdu
            if gca_param.plot_small
               blueMarker = '.';
               redMarker  = '^';
               blueMarkersize = 10; blueLineWidth = 2;
               redMarkersize  = 4;  redLineWidth  = 1;
            else
               blueMarker = '^';
               redMarker  = '^';
               blueMarkersize = 6;  blueLineWidth = 1;
               redMarkersize  = 6;  redLineWidth  = 1;
            end
         elseif ( allInd(iii,2)-allInd(iii,1)+1 ) == num_sim_norm
            if gca_param.plot_small
               blueMarker = '.';
               redMarker  = 'x';
               blueMarkersize = 10; blueLineWidth = 2;
               redMarkersize  = 4;  redLineWidth  = 1;
            else
               blueMarker = 'x';
               redMarker  = 'x';
               blueMarkersize = 6;  blueLineWidth = 1.5;
               redMarkersize  = 6;  redLineWidth  = 1.5;
            end
         else
            if gca_param.plot_small
               blueMarker = '.';
               redMarker  = 's';
               blueMarkersize = 10; blueLineWidth = 2;
               redMarkersize  = 4;  redLineWidth  = 1;
            else
               blueMarker = 's';
               redMarker  = 's';
               blueMarkersize = 6;  blueLineWidth = 1;
               redMarkersize  = 6;  redLineWidth  = 1;
            end
         end
         if iii > 1
            iii_offset = allInd(iii-1,2);
         else
            iii_offset = 0;
         end
         Iblue = iii_offset + find( allG(1,allInd(iii,1):allInd(iii,2))-g0 > 0 );
         Ired  = iii_offset + find( allG(1,allInd(iii,1):allInd(iii,2))-g0 < 0 );
         h1 = plot(allU(1,Iblue),allU(2,Iblue),blueMarker);
         set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
         h2 = plot(allU(1,Ired),allU(2,Ired),redMarker);
         set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
      end

      if n_pt_supp_clust > 0
         h_clust = plot(U_clust(1,:),U_clust(2,:),'ko','Markersize',8,'MarkerFaceColor','r');
      end
      if n_pt_supp_close > 0
         if flag_cluster_supp_close
            h_close0 = plot(U_close0(1,:),U_close0(2,:),'ks','Markersize',8,'MarkerFaceColor','g');
         end
         h_close = plot(U_close(1,:),U_close(2,:),'ko','Markersize',8,'MarkerFaceColor','g');
      end
      if n_pt_supp_switch > 0
         h_switch = plot(U_switch(1,:),U_switch(2,:),'ko','Markersize',8,'MarkerFaceColor','b');
      end

      if gca_param.auto_gca
         set(gca,'FontSize',14);
         axis equal
      else
         set(gca,'FontSize',14,'XLim',gca_param.XLim,'YLim',gca_param.YLim,'XTick',gca_param.XTick,'YTick',gca_param.YTick);
         axis square
      end
      grid on
      box on

      xlabel(probdata.name(1),'FontSize',14);
      ylabel(probdata.name(2),'FontSize',14);

      if analysisopt.flag_pause
         disp('Press a key'), pause
      else
         pause(1)
      end

   elseif nrv==3 & flag_plot23_all

      figure(20+Nb_step)
      clf
      if gca_param.plot_thresholds
         for istep = 0:(Nb_step-1)
            G0(istep+1) = svmdata(istep+1).g0;
         end
         G0(Nb_step+1) = g0;
         if n_pt_supp_clust > 0
            plot_svm3d(SVC,gca_param, ...
                       Usvm,n_pt_supp_clust,Usvm_in_margin,[],n_pt_supp_switch,Usvm_switch, ...
                       probdata,analysisopt,gfundata,Nb_step,G0);
         else
            plot_svm3d(SVC,gca_param, ...
                       Usvm,n_pt_supp_clust,[],Usvm_in_margin_sorted,n_pt_supp_switch,Usvm_switch, ...
                       probdata,analysisopt,gfundata,Nb_step,G0);
         end
      else
         plot_svm3d(SVC,gca_param, ...
                    Usvm,n_pt_supp_clust,Usvm_in_margin,Usvm_in_margin_sorted,n_pt_supp_switch,Usvm_switch);
      end

      if n_pt_supp_clust > 0
         for i_ball = 1:n_pt_supp_clust
            [X_ball,Y_ball,Z_ball] = sphere(20);
            X_ball = U_clust(1,i_ball)+1/5*X_ball;
            Y_ball = U_clust(2,i_ball)+1/5*Y_ball;
            Z_ball = U_clust(3,i_ball)+1/5*Z_ball;
            h_ball = surf(X_ball,Y_ball,Z_ball);
            set(h_ball,'FaceColor','r','EdgeColor','none');
         end
      end
      if n_pt_supp_close > 0
         for i_ball = 1:n_pt_supp_close
            [X_ball,Y_ball,Z_ball] = sphere(20);
            X_ball = U_close(1,i_ball)+1/5*X_ball;
            Y_ball = U_close(2,i_ball)+1/5*Y_ball;
            Z_ball = U_close(3,i_ball)+1/5*Z_ball;
            h_ball = surf(X_ball,Y_ball,Z_ball);
            set(h_ball,'FaceColor','g','EdgeColor','none');
         end
      end
      if n_pt_supp_switch > 0
         for i_ball = 1:n_pt_supp_switch
            [X_ball,Y_ball,Z_ball] = sphere(20);
            X_ball = U_switch(1,i_ball)+1/5*X_ball;
            Y_ball = U_switch(2,i_ball)+1/5*Y_ball;
            Z_ball = U_switch(3,i_ball)+1/5*Z_ball;
            h_ball = surf(X_ball,Y_ball,Z_ball);
            set(h_ball,'FaceColor','b','EdgeColor','none');
         end
      end
      
      if gca_param.auto_gca
         set(gca,'FontSize',14);
         axis equal
      else
         set(gca,'FontSize',14,'XLim',gca_param.XLim,'YLim',gca_param.YLim,'ZLim',gca_param.ZLim,'XTick',gca_param.XTick,'YTick',gca_param.YTick,'ZTick',gca_param.ZTick);
         axis square
      end
      grid on
      box on

      xlabel(probdata.name(1),'FontSize',14);
      ylabel(probdata.name(2),'FontSize',14);
      zlabel(probdata.name(3),'FontSize',14);

      if analysisopt.flag_pause
         disp('Press a key'), pause
      else
         pause(1)
      end
      
      
   end

   % Evaluation of ( n_pt_supp_clust + n_pt_supp_close + n_pt_supp_switch ) points (calls to the limit-state function)
   num_sim = ( n_pt_supp_clust + n_pt_supp_close + n_pt_supp_switch );

   k = 0;
   block_size = analysisopt.block_size;
   while k < num_sim

      block_size = min( block_size, num_sim-k );
      k = k + block_size;

      block_u = U(:,(k-block_size+1):k);

      block_x = u_to_x(block_u,probdata);

      [ block_g, dummy ] = gfun(lsf,block_x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);

      allU = [ allU block_u ];
      allG = [ allG block_g ];

   end

   allInd = [ allInd; Nfun+1 Nfun+num_sim ];
   Nfun = Nfun + num_sim;

   if Nb_step == 0

      I_train = 1:length(allG);

   else

      gridsel_option_svc = 0;
      num_sim = length(allG);
      k = 0;
      I_train = [];
      block_size = analysisopt.buf_size;
      svm_buf_size = analysisopt.svm_buf_size;
      while k < num_sim
         block_size = min( block_size, num_sim-k );
         k = k + block_size;
         block_u = allU(:,(k-block_size+1):k);
         block_g = eval_svm(block_u,svm_buf_size,SVC0,gridsel_option_svc);
         I_train = [ I_train find( block_g < 0 ) ];
      end

   end

   % Update of SVM classifier, based on the num_sim = ( n_pt_supp_clust + n_pt_supp_close + n_pt_supp_switch ) new training points
   if Nb_step == 0
      if kk < n_supp2
         if gridsel_svc_flag
            gridsel_option_svc = 1;
         else
            gridsel_option_svc = 0;
         end
         rbf_tab = rbf_tab0;
      elseif kk < n_supp3
         if gridsel_svc_flag
            gridsel_option_svc = 1;
            if (ind_SVC0+ind_SVC_best) > 1 & (ind_SVC0+ind_SVC_best) < length(rbf_tab0)
               rbf_tab = rbf_tab0([ind_SVC0+ind_SVC_best-1 ind_SVC0+ind_SVC_best ind_SVC0+ind_SVC_best+1]);
               ind_SVC0 = ind_SVC0+ind_SVC_best-2;
            elseif (ind_SVC0+ind_SVC_best) <= 1
               rbf_tab = rbf_tab0([1 2]);
               ind_SVC0 = 0;
            elseif (ind_SVC0+ind_SVC_best) >= length(rbf_tab0)
               rbf_tab = rbf_tab0([length(rbf_tab0)-1 length(rbf_tab0)]);
               ind_SVC0 = length(rbf_tab0)-2;
            end
         else
            gridsel_option_svc = 0;
         end
      else
         if gridsel_svc_flag
            gridsel_option_svc = 0;
            rbf_tab = rbf_tab0(ind_SVC0+ind_SVC_best);
         else
            gridsel_option_svc = 0;
         end
      end
   else
      if flag_var_rbf == 1
         if kk > n_supp2 & kk < n_supp3
            if gridsel_svc_flag
               gridsel_option_svc = 1;
               if (ind_SVC0+ind_SVC_best) > 1 & (ind_SVC0+ind_SVC_best) < length(rbf_tab0)
                  rbf_tab = rbf_tab0([ind_SVC0+ind_SVC_best-1 ind_SVC0+ind_SVC_best ind_SVC0+ind_SVC_best+1]);
                  ind_SVC0 = ind_SVC0+ind_SVC_best-2;
               elseif (ind_SVC0+ind_SVC_best) <= 1
                  rbf_tab = rbf_tab0([1 2]);
                  ind_SVC0 = 0;
               elseif (ind_SVC0+ind_SVC_best) >= length(rbf_tab0)
                  rbf_tab = rbf_tab0([length(rbf_tab0)-1 length(rbf_tab0)]);
                  ind_SVC0 = length(rbf_tab0)-2;
               end
            else
               gridsel_option_svc = 0;
            end
         else
            if gridsel_svc_flag
               gridsel_option_svc = 0;
               rbf_tab = rbf_tab0(ind_SVC0+ind_SVC_best);
            else
               gridsel_option_svc = 0;
            end
         end
      elseif flag_var_rbf == 0
               if gridsel_svc_flag
                  gridsel_option_svc = 0;
                  rbf_tab = rbf_tab0(ind_SVC0+ind_SVC_best);
               else
                  gridsel_option_svc = 0;
               end
      end
   end
   SVC = svm_init(gridsel_option_svc,rbf_tab);
   Y = sign( allG(1,I_train) - g0 );
   D = data(allU(:,I_train)',Y')
   if ( Nb_step > 0 )
      tic;
   end
   [ r, SVC ] = train(SVC,D);
   if ( Nb_step > 0 )
      training_time = toc;
      fprintf(1,'Training time: %f seconds\n', training_time);
   end
   if gridsel_option_svc == 1
      ind_SVC_best = SVC.best_index;
      gridsel_option_svc = 2;
      SVC = svm_init(gridsel_option_svc,rbf_tab,ind_SVC_best);
      [ r, SVC ] = train(SVC,D);
      rbf(kk) = rbf_tab(ind_SVC_best);
   else
      rbf(kk) = rbf_tab;
   end
   
   ssvm_misclassified

end % while kk < n_supp_stop

% Plots for nrv==2 or nrv==3
if nrv==2 & ( flag_plot23_all | flag_plot23_lastonly )
   
   figure(20+Nb_step)
   clf
   if gca_param.plot_thresholds
      for istep = 0:(Nb_step-1)
         G0(istep+1) = svmdata(istep+1).g0;
      end
      G0(Nb_step+1) = g0;
      plot_svm2d(SVC,gca_param, ...
                 probdata,analysisopt,gfundata,Nb_step,G0);
   else
      plot_svm2d(SVC,gca_param);
   end

   for iii = 1:size(allInd,1)
      if ( allInd(iii,2)-allInd(iii,1)+1 ) == num_sim_sdu
         if gca_param.plot_small
            blueMarker = '.';
            redMarker  = '^';
            blueMarkersize = 10; blueLineWidth = 2;
            redMarkersize  = 4;  redLineWidth  = 1;
         else
            blueMarker = '^';
            redMarker  = '^';
            blueMarkersize = 6;  blueLineWidth = 1;
            redMarkersize  = 6;  redLineWidth  = 1;
         end
      elseif ( allInd(iii,2)-allInd(iii,1)+1 ) == num_sim_norm
         if gca_param.plot_small
            blueMarker = '.';
            redMarker  = 'x';
            blueMarkersize = 10; blueLineWidth = 2;
            redMarkersize  = 4;  redLineWidth  = 1;
         else
            blueMarker = 'x';
            redMarker  = 'x';
            blueMarkersize = 6;  blueLineWidth = 1.5;
            redMarkersize  = 6;  redLineWidth  = 1.5;
         end
      else
         if gca_param.plot_small
            blueMarker = '.';
            redMarker  = 's';
            blueMarkersize = 10; blueLineWidth = 2;
            redMarkersize  = 4;  redLineWidth  = 1;
         else
            blueMarker = 's';
            redMarker  = 's';
            blueMarkersize = 6;  blueLineWidth = 1;
            redMarkersize  = 6;  redLineWidth  = 1;
         end
      end
      if iii > 1
         iii_offset = allInd(iii-1,2);
      else
         iii_offset = 0;
      end
      Iblue = iii_offset + find( allG(1,allInd(iii,1):allInd(iii,2))-g0 > 0 );
      Ired  = iii_offset + find( allG(1,allInd(iii,1):allInd(iii,2))-g0 < 0 );
      h1 = plot(allU(1,Iblue),allU(2,Iblue),blueMarker);
      set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
      h2 = plot(allU(1,Ired),allU(2,Ired),redMarker);
      set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
   end

   if gca_param.auto_gca
      set(gca,'FontSize',14);
      axis equal
   else
      set(gca,'FontSize',14,'XLim',gca_param.XLim,'YLim',gca_param.YLim,'XTick',gca_param.XTick,'YTick',gca_param.YTick);
      axis square
   end
   grid on
   box on

   xlabel(probdata.name(1),'FontSize',14);
   ylabel(probdata.name(2),'FontSize',14);

   if analysisopt.flag_pause
      disp('Press a key'), pause
   else
      pause(1)
   end
   
elseif nrv==3 & ( flag_plot23_all | flag_plot23_lastonly )

   figure(20+Nb_step)
   clf
   if gca_param.plot_thresholds
      for istep = 0:(Nb_step-1)
         G0(istep+1) = svmdata(istep+1).g0;
      end
      G0(Nb_step+1) = g0;
      plot_svm3d(SVC,gca_param, ...
                 probdata,analysisopt,gfundata,Nb_step,G0);
   else
      plot_svm3d(SVC,gca_param);
   end

   if gca_param.auto_gca
      set(gca,'FontSize',14);
      axis equal
   else
      set(gca,'FontSize',14,'XLim',gca_param.XLim,'YLim',gca_param.YLim,'ZLim',gca_param.ZLim,'XTick',gca_param.XTick,'YTick',gca_param.YTick,'ZTick',gca_param.ZTick);
      axis square
   end
   grid on
   box on

   xlabel(probdata.name(1),'FontSize',14);
   ylabel(probdata.name(2),'FontSize',14);
   zlabel(probdata.name(3),'FontSize',14);

   if analysisopt.flag_pause
      disp('Press a key'), pause
   else
      pause(1)
   end

end


if n_supp_stop <= n_supp1
   Nusvm = [ nusvm1 ];
elseif n_supp_stop <= n_supp2
   Nusvm = [ nusvm1 nusvm2 ];
elseif n_supp_stop <= n_supp3
   Nusvm = [ nusvm1 nusvm2 nusvm3 ];
else
   if nusvm4 == nusvm3
      Nusvm = [ nusvm1 nusvm2 nusvm3 ];
   else
      Nusvm = [ nusvm1 nusvm2 nusvm3 nusvm4 ];
   end
end

svm_buf_size = analysisopt.svm_buf_size;

if g0 > 0

   fprintf(1,'\nFinal estimate for pf and storage of SVM germs\n');
   
   for i = 1:length(Nusvm)

      if Nb_step > 0
         
         germset = i;

         % Load germs
         Subgerm0 = [];
         Subgermg0 = [];
         for ik = 1:Ikmax(germset)
            fprintf(1,[ 'Germ from file step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat\n' ]);
            eval([ 'load step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat' ]);
            Subgerm0 = [ Subgerm0 subgerm0 ];
            Subgermg0 = [ Subgermg0 subgermg0 ];
         end

         %  Generation of num_sim = Nusvm(i) points(N1, N2 or N3) using the lambda rj-mM algorithm based on the current SVM classifier (please refer to Bourinet et al. Paper)
         num_sim = Nusvm(i);
         ratio_amp_factor = 1; % lambda factor = 1 - Original Modified Metropolis Hasting algorithm (termed rj-mM algorithm in Bourinet et al. paper)
         genpopsubset

      end

      % Evaluation of num_sim points based on the current SVM classifier
      if Nb_step == 0
         num_sim = Nusvm(i);
      else
         num_sim = size_pop_subset;
      end

      block_size = analysisopt.buf_size;
      k = 0;
      ik = 0;
      nb_neg       = 0;
      n_def_U      = 0;
      n_inf_minus1 = 0;
      n_inf_plus1  = 0;
      n_in_margin  = 0;
      n_U          = 0;
      while k < num_sim

         block_size = min( block_size, num_sim-k );
         k = k + block_size;
         ik = ik + 1;

         if Nb_step == 0
            block_u = inv_norm_cdf(twister(nrv,block_size));
         else
            eval([ 'load step' num2str(Nb_step) '_pop_' num2str(ik) '.mat' ]);
         end
         
         block_g = eval_svm(block_u,svm_buf_size,SVC,gridsel_option_svc);

         % germ with exactly 10% of the whole population
         I_neg = find( block_g <= 0 );
         if ( nb_neg + length(I_neg ) ) > ( 0.1*Nusvm(i) )
            if ~isempty( I_neg(1:(0.1*Nusvm(i)-nb_neg)) )
               subgerm0  = block_u(:,I_neg(1:(0.1*Nusvm(i)-nb_neg)));
               subgermg0 = block_g(I_neg(1:(0.1*Nusvm(i)-nb_neg)));
               eval([ 'save step' num2str(Nb_step) '_subgerm' num2str(i) '_' num2str(ik) '.mat subgerm0 subgermg0' ]);
               Ikmax(i) = ik;
            end
         else
            subgerm0  = block_u(:,I_neg);
            subgermg0 = block_g(I_neg);
            eval([ 'save step' num2str(Nb_step) '_subgerm' num2str(i) '_' num2str(ik) '.mat subgerm0 subgermg0' ]);
            Ikmax(i) = ik;
         end
         nb_neg = nb_neg + length(I_neg);
         
         if i == length(Nusvm)
            n_def_U      = n_def_U + length(I_neg);                        % Number of failure points
            n_inf_minus1 = n_inf_minus1 + length(find( block_g <= -1 ));   % Number of points < -1
            n_inf_plus1  = n_inf_plus1 + length(find( block_g <= 1 ));     % Number of points < +1
            n_in_margin  = n_in_margin + length(find( 1 > abs(block_g) )); % Number of points in the margin
            n_U          = n_U + length( block_g );                        % Total number of points
         end

      end

   end
   
else

   fprintf(1,'\nFinal estimate for pf\n');
   
   for i = length(Nusvm):length(Nusvm)

      if Nb_step > 0
         
         germset = i;

         % Load germs
         Subgerm0 = [];
         Subgermg0 = [];
         for ik = 1:Ikmax(germset)
            fprintf(1,[ 'Germ from file step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat\n' ]);
            eval([ 'load step' num2str(Nb_step-1) '_subgerm' num2str(germset) '_' num2str(ik) '.mat' ]);
            Subgerm0 = [ Subgerm0 subgerm0 ];
            Subgermg0 = [ Subgermg0 subgermg0 ];
         end

         % Generation of num_sim = Nusvm(i) points using the modified Metropolis-Hastings algorithm based on the current SVM classifier
         num_sim = Nusvm(i);
         ratio_amp_factor = 1;
         genpopsubset

      end

      % Evaluation of num_sim points based on the current SVM classifier
      if Nb_step == 0
         num_sim = Nusvm(i);
      else
         num_sim = size_pop_subset;
      end

      block_size = analysisopt.buf_size;
      k = 0;
      ik = 0;
      nb_neg       = 0;
      n_def_U      = 0;
      n_inf_minus1 = 0;
      n_inf_plus1  = 0;
      n_in_margin  = 0;
      n_U          = 0;
      while k < num_sim

         block_size = min( block_size, num_sim-k );
         k          = k + block_size;
         ik         = ik + 1;

         if Nb_step == 0
            block_u = inv_norm_cdf(twister(nrv,block_size));
         else
            eval([ 'load step' num2str(Nb_step) '_pop_' num2str(ik) '.mat' ]);
         end
         
         block_g = eval_svm(block_u,svm_buf_size,SVC,gridsel_option_svc);

         I_neg = find( block_g <= 0 );
         if i == length(Nusvm)
            n_def_U      = n_def_U + length(I_neg);                        % Number of failure points
            n_inf_minus1 = n_inf_minus1 + length(find( block_g <= -1 ));   % Number of points < -1
            n_inf_plus1  = n_inf_plus1 + length(find( block_g <= 1 ));     % Number of points < +1
            n_in_margin  = n_in_margin + length(find( 1 > abs(block_g) )); % Number of points in the margin
            n_U          = n_U + length( block_g );                        % Total number of points
         end
                  
      end
      
   end
   
end

fprintf(1,'\nStep%d - Final evaluation                 :\n', Nb_step);
fprintf(1,'\nNumber of failure points             : %d\n', n_def_U);
fprintf(1,'Total number of learning points      : %d\n', n_U)
fprintf(1,'Pf estimate                          : %e\n', n_def_U/n_U);
fprintf(1,'Pf < [ -1 , 0 , +1 ]                 : [ %e , %e , %e ]\n', n_inf_minus1/n_U,n_def_U/n_U,n_inf_plus1/n_U);
svmpf(1,kk+1) = n_def_U/n_U;
svmpf(2,kk+1) = n_inf_minus1/n_U;
svmpf(3,kk+1) = n_inf_plus1/n_U;
fprintf(1,'Number of points in the margin       : %d\n',n_in_margin);
fprintf(1,'Ratio of points in the margin / total: %f\n\n',n_in_margin/n_U);

svmdata(Nb_step+1).n_supp1          = n_supp1;
svmdata(Nb_step+1).n_supp2          = n_supp2;
svmdata(Nb_step+1).n_supp3          = n_supp3;
svmdata(Nb_step+1).n_supp_stop      = n_supp_stop;
svmdata(Nb_step+1).g0               = g0;
svmdata(Nb_step+1).bestIter         = bestIter;
svmdata(Nb_step+1).SVC              = SVC;
svmdata(Nb_step+1).rbf              = rbf;
svmdata(Nb_step+1).svmpf            = svmpf;
svmdata(Nb_step+1).svmmargin        = svmmargin;
svmdata(Nb_step+1).svmswitch        = svmswitch;
svmdata(Nb_step+1).svmmisclassified = svmmisclassified;