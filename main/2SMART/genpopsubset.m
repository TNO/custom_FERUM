% Generate num_sim points using the lambda rj-mM algorithm based on the current SVM classifier (please refer to Bourinet et al. paper for details)
% For lambda = 1, lambda rj-mM algorithm is merely the original modified Metropolis-Hastings algorithm 

buf_size         = analysisopt.buf_size;

subset           = [];
subsetg          = [];
size_pop_subset  = 0;
threshold_subset = 0;
Nb_generation    = 0;
ik               = 0;

while size_pop_subset < num_sim
   
   Nb_generation = Nb_generation + 1;
   fprintf(1,'Generation #%d\n',Nb_generation)
   
   ssvm_MCMC
   
   % Initializations
   num_sim_generation = size(Subtemp,2);
   k                  = 0;
   block_size         = analysisopt.buf_size;
   svm_buf_size       = analysisopt.svm_buf_size;
   Subtempg           = zeros(1,num_sim_generation);
   
   while k < num_sim_generation
      block_size = min( block_size, num_sim_generation-k );
      k = k + block_size;
      block_u = Subtemp(:,(k-block_size+1):k);
      block_g = eval_svm(block_u,svm_buf_size,SVC0,0);
      Subtempg((k-block_size+1):k) = block_g;
   end % while k < num_sim_generation
   
   ind = find( Subtempg > 0 );

   Subtemp(:,ind) = Subgerm0(:,ind);
   Subtempg(ind)  = Subgermg0(ind);
      
   if Nb_generation > 10
      
      subset  = [ subset Subtemp ];
      subsetg = [ subsetg Subtempg ];
      
      size_subset = length(subsetg);
  
      ii = 0;
      while size_subset >= buf_size

         ii = ii + 1;

         if ( size_pop_subset + buf_size ) < num_sim
            
            block_u =  subset( : , ((ii-1)*buf_size+1) : ((ii-1)*buf_size+buf_size) );
            block_g = subsetg( 1 , ((ii-1)*buf_size+1) : ((ii-1)*buf_size+buf_size) );
            size_subset     = size_subset     - buf_size;
            size_pop_subset = size_pop_subset + buf_size;
            subset( : ,  ((ii-1)*buf_size+1) : ((ii-1)*buf_size+buf_size) ) = nan;
            subsetg( 1 , ((ii-1)*buf_size+1) : ((ii-1)*buf_size+buf_size) ) = nan;

         else
            
            block_u =  subset( : , ((ii-1)*buf_size+1) : ((ii-1)*buf_size+(num_sim-size_pop_subset)) );
            block_g = subsetg( 1 , ((ii-1)*buf_size+1) : ((ii-1)*buf_size+(num_sim-size_pop_subset)) );
            size_subset     = size_subset     - ( num_sim - size_pop_subset );
            size_pop_subset = size_pop_subset + ( num_sim - size_pop_subset );
            subset( : ,  ((ii-1)*buf_size+1) : ((ii-1)*buf_size+(num_sim-size_pop_subset)) ) = nan;
            subsetg( 1 , ((ii-1)*buf_size+1) : ((ii-1)*buf_size+(num_sim-size_pop_subset)) ) = nan;
            
         end
            
         ik = ik + 1;
         fprintf(1,[ 'save step' num2str(Nb_step) '_pop_' num2str(ik) '.mat block_u block_g\n' ]);
         eval([ 'save step' num2str(Nb_step) '_pop_' num2str(ik) '.mat block_u block_g' ]);

         if size_pop_subset == num_sim, break, end
         
      end
      
      if ( size_subset < buf_size ) & ( size_subset > num_sim )

         block_u =  subset( : , 1 : num_sim );
         block_g = subsetg( 1 , 1 : num_sim );
         size_subset     = size_subset     - num_sim;
         size_pop_subset = size_pop_subset + num_sim;
         subset( :  , 1 : num_sim ) = nan;
         subsetg( 1 , 1 : num_sim ) = nan;

         ik = ik + 1;
         fprintf(1,[ 'save step' num2str(Nb_step) '_pop_' num2str(ik) '.mat block_u block_g\n' ]);
         eval([ 'save step' num2str(Nb_step) '_pop_' num2str(ik) '.mat block_u block_g' ]);

      end
      
      Jnan = find(isnan(subsetg));
      if ~isempty(Jnan)
         subset(:,Jnan) = [];
         subsetg(Jnan)  = [];
      end
     
   end % if Nb_generation > 10
   
   Subgerm0  = Subtemp;
   Subgermg0 = Subtempg;
   
end % while size_pop_subset < num_sim