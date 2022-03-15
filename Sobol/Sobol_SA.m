function [ svrdata, sobolresults, probdata ] = Sobol_SA(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun
nfun = 0;

% Options settings

marg           = probdata.marg;
R              = probdata.correlation;
transf_type    = probdata.transf_type;
flag_sens      = probdata.flag_sens;
nrv            = size(marg,1);

NbCal          = analysisopt.NbCal;

sampling       = analysisopt.sampling;
block_size     = analysisopt.block_size;
rand_generator = analysisopt.rand_generator;
num_sim        = analysisopt.num_sim;
first_indices  = analysisopt.first_indices;
total_indices  = analysisopt.total_indices;
all_indices    = analysisopt.all_indices;

if isfield(analysisopt,'SVR')
   optSVR     = analysisopt.SVR;
   SVRbasis   = analysisopt.SVRbasis;
   SVR_Nbasis = analysisopt.SVR_Nbasis;
else
   optSVR = 'no';
end

if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end

if isfield(gfundata(lsf),'ng')
   ng = gfundata(lsf).ng;
else
   ng = 1;
end

% Initialization

k = 0;
i = 0;

S  = [];
ST = [];

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
      Lo = (chol(Ro))'; probdata.Lo = Lo;
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


% Evaluation of indices on the model directly or on the SVR

if strcmp(optSVR ,'yes') % Build a SVR model

   gfundata = train_SVR(lsf,probdata,analysisopt,gfundata,femodel,randomfield);
   
   sobolresults.nfun = nfun;
      
   for ig = 1:ng

      % SVR parameters

      svrdata(ig).hyperparameters = gfundata(lsf).SVR(ig).hyperparameters;
      svrdata(ig).error           = gfundata(lsf).SVR(ig).error;

   end

   initial_block_size = analysisopt.block_size;
   analysisopt.block_size = analysisopt.buf_size;   % Number of calls to the g SVR-surrogate to be sent simultaneously

   if NbCal == 1

      [ retour, retourT, retourTout ] = Sobol_evaluation(1,probdata,analysisopt,gfundata,femodel,randomfield);

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 0
         sobolresults.First = retour;
      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1
         sobolresults.FirstSecond = retour;
      end

      if analysisopt.total_indices == 1
         sobolresults.Total = retourT;
      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
         sobolresults.Total = retourTout;
      end

      sobolresults.Total = retourT;
      sobolresults.All   = retourTout;

   else

      for k = 1:NbCal

         [ retour, retourT, retourTout ] = Sobol_evaluation(1,probdata,analysisopt,gfundata,femodel,randomfield);

         if analysisopt.first_indices == 1 & analysisopt.all_indices == 0
            First(:,k,1:ng) = retour(:,1+[1:ng]);
         end
         
         if analysisopt.first_indices == 1 & analysisopt.all_indices == 1
            FirstSecond(:,k,1:ng) = retour(:,1+[1:ng]);
         end
         
         if analysisopt.total_indices == 1
            Total(:,k,1:ng) = retourT(:,1+[1:ng]);
         end
         
         if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
            All(:,k,1:ng) = retourTout(:,1+[1:ng]);
         end

      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 0
         for ig = 1:ng
            first(:,:,ig) = First(:,:,ig)';
         end
         ResFirst(:,1) = retour(:,1);
         for ig = 1:ng
            ResFirst(:,1+ig)    = mean(first(:,:,ig))';
            ResFirst(:,1+ng+ig) = sqrt(var(first(:,:,ig)))';
         end
         sobolresults.First = ResFirst;
      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1
         for ig = 1:ng
            firstsecond(:,:,ig) = FirstSecond(:,:,ig)';
         end
         ResFirstSecond(:,1) = retour(:,1);
         for ig = 1:ng
             ResFirstSecond(:,1+ig)    = mean(firstsecond(:,:,ig))';
             ResFirstSecond(:,1+ng+ig) = sqrt(var(firstsecond(:,:,ig)))';
         end
         sobolresults.FirstSecond = ResFirstSecond;
      end

      if analysisopt.total_indices == 1
         for ig = 1:ng
            total(:,:,ig) = Total(:,:,ig)';
         end
         ResTotal(:,1) = retourT(:,1);
         for ig = 1:ng
             ResTotal(:,1+ig)    = mean(total(:,:,ig))';
             ResTotal(:,1+ng+ig) = sqrt(var(total(:,:,ig)))';
         end
         sobolresults.Total = ResTotal;
      end
      
      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
         for ig = 1:ng
            all(:,:,ig) = All(:,:,ig)';
         end
         ResAl(:,1) = retourTout(:,1);
         for ig = 1:ng
             ResAll(:,1+ig)    = mean(all(:,:,ig))';
             ResAll(:,1+ng+ig) = sqrt(var(all(:,:,ig)))';
         end
         sobolresults.All = ResAll;
      end

   end

   analysisopt.block_size = initial_block_size;

else  % Sobol' indices evaluation on the original model

   svrdata.exist = 'no';

   if NbCal == 1

      [ retour, retourT, retourTout ] = Sobol_evaluation(1,probdata,analysisopt,gfundata,femodel,randomfield);

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 0
         sobolresults.First = retour;
      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1
         sobolresults.FirstSecond = retour;
      end
      
      if analysisopt.total_indices == 1
         sobolresults.Total = retourT;
      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
         sobolresults.Total = retourTout;
      end

      sobolresults.Total = retourT;
      sobolresults.All   = retourTout;
      
      sobolresults.nfun  = nfun;


   else
      
      for k = 1:NbCal
         
         [ retour, retourT, retourTout ] = Sobol_evaluation(1,probdata,analysisopt,gfundata,femodel,randomfield);

         if analysisopt.first_indices == 1 & analysisopt.all_indices == 0
            First(:,k,1:ng) = retour(:,1+[1:ng]);
         end
         
         if analysisopt.first_indices == 1 & analysisopt.all_indices == 1
            FirstSecond(:,k,1:ng) = retour(:,1+[1:ng]);
         end
         
         if analysisopt.total_indices == 1
            Total(:,k,1:ng) = retourT(:,1+[1:ng]);
         end
         
         if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
            All(:,k,1:ng) = retourTout(:,1+[1:ng]);
         end
         
      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 0
         for ig = 1:ng
            first(:,:,ig) = First(:,:,ig)';
         end
         ResFirst(:,1) = retour(:,1);
         for ig = 1:ng
            ResFirst(:,1+ig)    = mean(first(:,:,ig))';
            ResFirst(:,1+ng+ig) = sqrt(var(first(:,:,ig)))';
         end
         sobolresults.First = ResFirst;
      end

      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1
         for ig = 1:ng
            firstsecond(:,:,ig) = FirstSecond(:,:,ig)';
         end
         ResFirstSecond(:,1) = retour(:,1);
         for ig = 1:ng
             ResFirstSecond(:,1+ig)    = mean(firstsecond(:,:,ig))';
             ResFirstSecond(:,1+ng+ig) = sqrt(var(firstsecond(:,:,ig)))';
         end
         sobolresults.FirstSecond = ResFirstSecond;
      end

      if analysisopt.total_indices == 1
         for ig = 1:ng
            total(:,:,ig) = Total(:,:,ig)';
         end
         ResTotal(:,1) = retourT(:,1);
         for ig = 1:ng
             ResTotal(:,1+ig)    = mean(total(:,:,ig))';
             ResTotal(:,1+ng+ig) = sqrt(var(total(:,:,ig)))';
         end
         sobolresults.Total = ResTotal;
      end
      
      if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
         for ig = 1:ng
            all(:,:,ig) = All(:,:,ig)';
         end
         ResAl(:,1) = retourTout(:,1);
         for ig = 1:ng
             ResAll(:,1+ig)    = mean(all(:,:,ig))';
             ResAll(:,1+ng+ig) = sqrt(var(all(:,:,ig)))';
         end
         sobolresults.All = ResAll;
      end
      
      sobolresults.nfun = nfun;
      
   end
   
end