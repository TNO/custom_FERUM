function data = subset_subset_step(data,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield)

% Perform a single subset step

global nfun

nrv = size(probdata.marg,1);

rand_generator = analysisopt.rand_generator;
if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end
flag_plot      = analysisopt.flag_plot;
flag_plot_gen  = analysisopt.flag_plot_gen;
 
width    = data.width;
Indgerm  = data.Indgerm(end,1:length(find(data.Indgerm(end,:))));
subgermU = data.U(:,Indgerm);
subgermG = data.G(Indgerm);  

subsetU   = [];
subsetG  = [];

Nb_generation = 0;
MHrate = [];
mMHrate = [];

while size(subsetU,2) < num_sim
  
   Nb_generation = Nb_generation + 1;

   [ subtempu, I_eval, newG, ratio ] = subset_MCMC(subgermU,width,rand_generator);
   MHrate = [ MHrate ; 100*mean(ratio,2) 100*sum(newG,2)/length(newG) ];

   % Initializations
   
   block_size = analysisopt.block_size;
   
   k = 0;
   percent_done = 0;

   num_subtemp_eval = length(I_eval);
   subtempu_eval = subtempu(:,I_eval);
   subtempg = subgermG;
   subtempg_eval = subtempg(1,I_eval);
      
   while k < num_subtemp_eval
   
      block_size = min( block_size, num_subtemp_eval-k );
      k = k + block_size;
      
      allu = subtempu_eval(:,(k-block_size+1):k);

      % Transform into original space
      allx = u_to_x(allu,probdata);

      % Evaluate limit-state function
      [ allG, dummy ] = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
      
      subtempg_eval((k-block_size+1):k) = allG;

      if floor( k/num_subtemp_eval * 20 ) > percent_done
         percent_done = floor( k/num_subtemp_eval * 20 );
        if echo_flag
            fprintf(1,'Subset step #%d - Generation #%d - %d%% complete\n',data.Nb_step,Nb_generation,percent_done*5)
        end
      end
      
   end
   
   subtempg(1,I_eval) = subtempg_eval;
   
   data.Neval = data.Neval + num_subtemp_eval;
   data.N     = data.N     + size(subtempu,2);
   
   % Selection des points dont la valeur de l'état-limite est supérieure au seuil en cours
   ind = find( subtempg > data.y(end) );

   if nrv == 2 & flag_plot == 1 & flag_plot_gen == 1
      data.Usubtemp = subtempu;
      data.ind = ind;
      data.Nb_generation = Nb_generation;
      data = subset_subset_graph(data,2);
   end
  
   mMHrate = [ mMHrate ; 100*(length(subtempg)-length(ind))/length(subtempg) ];

   % On remplace les points selectionnés par les anciens points
   subtempu(:,ind) = subgermU(:,ind);
   subtempg(ind) = subgermG(ind);

   subsetU = [ subsetU subtempu ]; 
   subsetG = [ subsetG subtempg ];

   subgermG = subtempg;
   subgermU = subtempu;
  
end

data.U       = [ data.U subsetU ];
data.G       = [ data.G subsetG ];
data.Indices = [ data.Indices; data.Indices(end,2)+1 data.Indices(end,2)+size(subsetU,2) ];
data.AccRate = [ data.AccRate [ Nb_generation ; mean(MHrate(:,1)) ; mean(MHrate(:,2)) ; mean(mMHrate) ] ];
