%     Finite Element Reliability Using Matlab, FERUM, Version 4.0, 2009 
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     A copy of the GNU General Public License is found in the file
%     <gpl-3.0.txt> following this collection of FERUM program files.
%     This license can be also found at: <http://www.gnu.org/licenses/>.    
%
%     For more information on FERUM, visit: <http://www.ifma.fr/FERUM/>


if exist('analysisopt')

   if isfield(analysisopt,'echo_flag')
      echo_flag = analysisopt.echo_flag;
   else
      echo_flag = 1;
   end
   
end


if echo_flag

   clc, format short e

   disp('   ________________________________________________________________________');
   disp('  | Welcome to FERUM Version 4.0 (Finite Element Reliability Using Matlab) |');
   disp('  | For more information, visit: http://www.ifma.fr/FERUM/                 |');
   disp('  | Note: All the analysis options below assumes that necessary data       |');
   disp('  |       are available in the current Matlab workspace.                   |');
   disp('   ________________________________________________________________________');
   disp(' ');
   disp('   0: Exit');
   disp('   1: Help');
   disp('  10: FORM Analysis');
   disp('  11: FORM Analysis - Multiple design point');
   disp('  12: SORM Analysis - Curvature Fitting method with computation of the Hessian');
   disp('  13: SORM Analysis - Point Fitting method');
   disp('  20: Distribution analysis');
   disp('  21: Importance Sampling / Monte Carlo Simulation Analysis');
   disp('  22: Directional Simulation Analysis');
   disp('  23: Subset Simulation Analysis');
   disp('  33: 2SMART Analysis');
   disp('  40: RBDO - Smooth Nested Bi-Level Approach (using the gradient of Pf w.r.t. thetag)');
   disp('  50: Sobol'' Global Sensitivity Analysis');
   disp(' ');
   analysistype = input('  CHOOSE OPTION FROM THE LIST ABOVE: ');
   analysisopt.analysistype = analysistype;

end


if analysisopt.analysistype > 2
   
   if ~isfield(analysisopt,'already_updated')
      % This function updates probdata and gfundata before any analysis (must be run only once)
      if echo_flag
         [probdata,gfundata,analysisopt] = update_data(1,probdata,analysisopt,gfundata,femodel), dummy=input('Hit a key to continue\n');
      else
         [probdata,gfundata,analysisopt] = update_data(1,probdata,analysisopt,gfundata,femodel);
      end
   end
   
   % This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
   probdata.marg  = distribution_parameter(probdata.marg);
   
end


tic


switch analysisopt.analysistype
   
   
   
   
   case 0 % ---- EXIT ----------------------------------------------------------------------------
      
      disp(' ');
      disp('  Bye, bye.');
      disp(' ');
   
   
      
   
   case 1 % ---- HELP ----------------------------------------------------------------------------

      clc
      disp(' ');
      disp('  FERUM HELP');
      disp(' ');
      disp('  If you are new to FERUM, you are recommended to visit the web page:');
      disp('  http://www.ifma.fr/FERUM');
      disp(' ');
      disp('  To run FERUM, do the following:');
      disp('  1. Specify necessary parameters in your current Matlab workspace. ');
      disp('     (The format of the input can be found in the inputfile_template.m file.');
      disp('     If you want to try one of the provided example inputfiles,');
      disp('     simply read the file into your workspace by writing the file name');
      disp('     without .m extension and press enter. ');
      disp('  2. Start FERUM (the shell program) by issuing the command ferum.');
      disp('  3. Choose the alternative that fits your purpose.');
      disp('  4. Main results as displayed on the screen and more detailed results are usually');
      disp('     stored in a specific data structure (e.g. formresults for a FORM analysis).');
      disp(' ');
   
   
      
   
   case 10 % ---- FORM ----------------------------------------------------------------------------
   
      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('FORM analysis is running, please wait... (Ctrl+C breaks)')
      end % End if echo_flag
   
      % Run the analysis
      [ formresults, probdata ] = form(1,probdata,analysisopt,gfundata,femodel,randomfield);

      % Store formresults in analysisopt structure for further potential IS simulations
      analysisopt.formresults = formresults;

      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING FORM RELIABILITY ANALYSIS' ])
         disp([' '])
         disp('Number of iterations:'), disp(formresults.iter)
         disp('Reliability index beta:'), disp(formresults.beta)
         disp('Failure probability:'), disp(formresults.pf1)
         disp('Number of calls to the limit-state function: '), disp(formresults.nfun)
         disp(['..............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   formresults.iter          = Number of iterations')
         disp('   formresults.beta          = Reliability index beta from FORM analysis')
         disp('   formresults.pf1           = Failure probability pf1')
         disp('   formresults.dsptu         = Design point u_star')
         disp('   formresults.alpha         = Alpha vector')
         disp('   formresults.dsptx         = Design point in original space')
         disp('   formresults.imptg         = Importance vector gamma')
         disp('   formresults.gfcn          = Recorded values of the limit-state function')
         disp('   formresults.stpsz         = Recorded step size values')
         disp('   formresults.e1            = Recorded values of g-criterion')
         disp('   formresults.e2            = Recorded values og u-criterion')
         disp('   formresults.Recorded_u    = Recorded values of u (optional)')
         disp('   formresults.Recorded_x    = Recorded values og x (optional)')
         if gfundata(1).flag_sens
            if isfield(gfundata(1),'thetag')
               disp('   formresults.dbeta_dthetag = Sensisitivities of beta index w.r.t. limit-state parameters')
               disp('   formresults.dpf_dthetag   = Sensisitivities of probability of failure w.r.t. limit-state parameters')
            end
         end
         switch probdata.flag_sens
            case 0
            otherwise
               disp('   formresults.dbeta_dthetaf = Sensitivities of beta index w.r.t. to distribution parameters')
               disp('   formresults.dpf_dthetaf   = Sensitivities of probability of failure w.r.t. to distribution parameters')
               disp('   formresults.dbeta_drho    = Sensitivities of beta index w.r.t. to correlation coefficients')
               disp('   formresults.dpf_drho      = Sensitivities of probability of failure w.r.t. to correlation coefficients')
               disp('   formresults.delta         = "Normalized" sensitivities of beta index w.r.t. to mean')
               disp('   formresults.eta           = "Normalized" sensitivities of beta index w.r.t. standard deviation')
         end
         disp('   formresults.nfun          = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
   
   
      
      
   case 11 % ---- Multiple design point FORM ------------------------------------------------------

      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('Multiple design point FORM analysis is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end % End if echo_flag

      % Run the analysis
      [ ALLformresults, probdata, analysisopt, gfundata ] = form_multiple_dspts(1,probdata,analysisopt,gfundata,femodel,randomfield);




   case 12 % ---- SORM - Curvature fitting method (Hessian) ---------------------------------------
      
      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('SORM analysis (Curvature Fitting method based on the computation of the hessian) is running,')
         disp('please wait... (Ctrl+C breaks)')
      end
   
      sormcfhresults = sorm_cfh(1,probdata,analysisopt,gfundata,femodel,randomfield);

      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SORM RELIABILITY ANALYSIS' ])
         disp([' '])
         disp('Reliability index beta (improved Breitung)'), disp(sormcfhresults.betag_breitung_m)
         disp('Failure probability:'), disp(sormcfhresults.pf2_breitung_m)
         disp('Number of calls to the limit-state function: '), disp(sormcfhresults.nfun)
         disp(['..............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   sormcfhresults.betag_breitung   = Reliability index beta from Breitung expression')
         disp('   sormcfhresults.pf2_breitung     = Failure probability from Breitung expression')
         disp('   sormcfhresults.betag_breitung_m = Reliability index beta from improved Breitung expression')
         disp('   sormcfhresults.pf2_breitung_m   = Failure probability from improved Breitung expression')
         disp('   sormcfhresults.kappa            = Curvatures at design point')
         disp('   sormcfhresults.nfun             = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
   
   
      
      
   case 13 % ---- SORM - Point fitting method -----------------------------------------------------

      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('SORM analysis (Point Fitting Method) is running, please wait... (Ctrl+C breaks)')
      end

      sormpfresults = sorm_pf(1,probdata,analysisopt,gfundata,femodel,randomfield);

      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SORM RELIABILITY ANALYSIS' ])
         disp([' '])
         disp('Reliability index beta (improved Breitung)'), disp(sormpfresults.betag_breitung_m)
         disp('Failure probability:'), disp(sormpfresults.pf2_breitung_m)
         disp('Number of calls to the limit-state function: '), disp(sormpfresults.nfun)
         disp(['..............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   sormpfresults.betag_breitung_m  = Reliability index beta from improved Breitung expression')
         disp('   sormpfresults.pf2_breitung_m    = Failure probability from improved Breitung expression')
         disp('   sormpfresults.uprimei_minus     = u''i- coordinates')
         disp('   sormpfresults.uprimen_minus     = u''n- ordinates')
         disp('   sormpfresults.G_minus           = G(u''-)')
         disp('   sormpfresults.ai_minus          = ai- coefficients')
         disp('   sormpfresults.uprimei_plus      = u''i+ coordinates')
         disp('   sormpfresults.uprimen_plus      = u''n+ ordinates')
         disp('   sormpfresults.G_plus            = G(u''+)')
         disp('   sormpfresults.ai_plus           = ai+ coefficients')
         disp('   sormpfresults.nfun              = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag

      
      
      
   case 20 % ---- DISTRIBUTION ANALYSIS -----------------------------------------------------------

      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('DISTRIBUTION analysis is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end

      % Run simulation analysis
      [ distributionresults, probdata ] = distribution_analysis(1,probdata,analysisopt,gfundata,femodel,randomfield);

      if echo_flag

         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING DISTRIBUTION ANALYSIS' ])
         disp([' '])
         disp('Number of simulations: '), disp(distributionresults.nfun)
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   distributionresults.Xmean  = Estimator of the means of simulated X variables')
         disp('   distributionresults.Xstdv  = Unbiased estimator of the variances of simulated X variables')
         disp('   distributionresults.Gmean  = Estimator of the mean of the simulated g outpout')
         disp('   distributionresults.Gstdv  = Unbiased estimator of the variance of the simulated g outpout')
         disp('   distributionresults.ALLX   = Simulated X variables')
         disp('   distributionresults.ALLG   = Simulated g outpout')
         disp('   distributionresults.nfun   = Number of simulations / calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])

      end % End if echo_flag
  
      
      
      
   case 21 % ---- CRUDE MONTE CARLO / IMPORTANCE SAMPLING SIMULATIONS -----------------------------

      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('CRUDE MONTE CARLO SIMULATION / IMPORTANCE SAMPLING is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end
   
      % Run simulation analysis
      if isfield(analysisopt,'lowRAM')
         lowRAM = analysisopt.lowRAM;
      else
         lowRAM = 0;
      end

      switch lowRAM
         case 0
            [ simulationresults, probdata ] = simulation_single_dspt(1,probdata,analysisopt,gfundata,femodel,randomfield);
         case 1
            [ simulationresults, probdata ] = simulation_single_dspt_lowRAM(1,probdata,analysisopt,gfundata,femodel,randomfield);
      end

      if echo_flag

         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING CRUDE MONTE CARLO SIMULATION / IMPORTANCE SAMPLING' ])
         disp([' '])
         disp('Probability of failure: '), disp(simulationresults.pf)
         disp('Coefficient of variation of the failure probability:'), disp(simulationresults.cov_pf)
         disp('Reliability index beta: '), disp(simulationresults.beta)
         disp('Number of simulations: '), disp(simulationresults.nfun)
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   simulationresults.pf     = Failure probability from Monte Carlo simulation')
         disp('   simulationresults.cov_pf = Coefficient of variation of the failure probability')
         disp('   simulationresults.beta   = Generalized reliability index beta from this simulation')
         disp('   simulationresults.nfun   = Number of simulations / calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])

      end % End if echo_flag




   case 22 % ---- DIRECTIONAL SIMULATIONS ---------------------------------------------------------
   
      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('DIRECTIONAL SIMULATION is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end
   
      % Run simulation analysis
      [ dirsimulationresults, probdata ] = dir_simulation(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING DIRECTIONAL SIMULATION' ])
         disp([' '])
         disp('Probability of failure: '), disp(dirsimulationresults.pf)
         if isfield(dirsimulationresults,'cov_pf')
               disp('Coefficient of variation of the failure probability:'), disp(dirsimulationresults.cov_pf)
         end
         disp('Reliability index beta: '), disp(dirsimulationresults.beta)
         disp('Number of directions used for simulations: '), disp(dirsimulationresults.ndir)
         disp('Number of calls to the limit-state function: '), disp(dirsimulationresults.nfun)
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   dirsimulationresults.pf     = Failure probability from directional Monte Carlo simulation')
         if isfield(dirsimulationresults,'cov_pf')
               disp('   dirsimulationresults.cov_pf = Coefficient of variation of the failure probability')
         end
         disp('   dirsimulationresults.beta   = Generalized reliability index beta from this simulation')
         disp('   dirsimulationresults.Alla   = Unit vectors of directions used for simultations')
         disp('                                 ( only when analysisopt.keep_a = 1 ) ')
         disp('   dirsimulationresults.Allr   = r values for which g(r) = 0 , r = rho if there is no intersection')
         disp('                                 ( only when analysisopt.keep_r = 1 )')
         disp('   dirsimulationresults.ndir   = Number of directions used for simulations')
         disp('   dirsimulationresults.nfun   = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag

      
      
      
   case 23 % ---- SUBSET SIMULATIONS ---------------------------------------------------------------

      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('SUBSET SIMULATION is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end

      % Run simulation analysis
      [ subsetsimulationresults, probdata ] = subset_simulation(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING SUBSET SIMULATION' ])
         disp([' '])
         disp('Probability of failure: '), disp(subsetsimulationresults.pf)
         if isfield(subsetsimulationresults,'cov_pf')
            disp('Coefficient of variation of failure probability (lower bound):'), disp(subsetsimulationresults.cov_pf)
         end
         disp('Reliability index beta: '), disp(subsetsimulationresults.beta)
         disp('Number of calls to the limit-state function: '), disp(subsetsimulationresults.nfun)
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   subsetsimulationresults.pf         = Failure probability from subset simulations')
         if isfield(subsetsimulationresults,'cov_pf')
            disp('   subsetsimulationresults.cov_pf     = Coefficient of variation for the failure probability (lower bound)')
         end
         disp('   subsetsimulationresults.beta       = Generalized reliability index beta from this simulation')
         disp('   subsetsimulationresults.SubsetData = Subset data structure')
         disp('   subsetsimulationresults.nfun       = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
      
      
      
      
   case 33 % ---- 2SMART ---------------------------------------------------------------------------
   
      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('2SMART analysis is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end
   
      % Run simulation analysis
      [ ssvmsimulationresults, probdata ] = ssvm_simulation(1,probdata,analysisopt,gfundata,femodel,randomfield);
   
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING 2SMART ANALYSIS' ])
         disp([' '])
         disp('Probability of failure: '), disp(ssvmsimulationresults.pf)
         disp('Reliability index beta: '), disp(ssvmsimulationresults.beta)
         disp('Number of calls to the limit-state function: '), disp(ssvmsimulationresults.nfun)
         disp(['...............................................................................................'])
         disp('The following parameters are now available in your current workspace:')
         disp('   ssvmsimulationresults.pf         = Failure probability from 2SMART analysis')
         disp('   ssvmsimulationresults.beta       = Generalized reliability index beta from this simulation')
         disp('   ssvmsimulationresults.SubsetData = Subset data structure')
         disp('   ssvmsimulationresults.nfun       = Number of calls to the limit-state function')
         disp(['..............................................................................................'])
         disp([' '])

      end % End if echo_flag




   case 40 % ---- N2-LA RBDO ANALYSIS based on gradient of beta w.r.t. design variables -----------
   
      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('N2LA analysis is running, please wait... (Ctrl+C breaks)')
         disp([' '])
      end
      
      % Run N2LA RBDO analysis
      [ N2LAresults, probdata, rbdo_parameters] = N2LA(1,rbdo_parameters,rbdo_fundata,probdata,analysisopt,gfundata,femodel,randomfield);
      
      if echo_flag
         
         % Display results
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         disp(['RESULTS FROM RUNNING N2-LA RBDO ANALYSIS' ])
         disp([' '])
         disp('Number of iterations:');disp(N2LAresults.number_of_iterations)
         disp([' '])
         disp('Number of calls to the limit-state function:');disp(N2LAresults.nfun)
         disp([' '])
         disp('Normalized optimal design point (thetag):');
         thetagname = gfundata(1).thetagname;
         nthetag = length(thetagname);
         for k=1:nthetag
            fprintf('     %s = %d\n',thetagname{k},N2LAresults.normalized_optimized_thetag(k,end))
         end
         disp([' '])
         disp('Optimal design point (thetag):');
         for k=1:nthetag
            fprintf('     %s = %d\n',thetagname{k},N2LAresults.optimized_thetag(k,end))
         end
         disp([' '])
         disp('Normalized optimized cost:');disp(N2LAresults.normalized_optimized_cost(end))
         disp('Optimized cost:');disp(N2LAresults.optimized_cost(end))
         disp([' '])
         disp('Reliability index:');disp(N2LAresults.optimized_beta(end))
         disp([' '])
         disp(['..............................................................................................'])
         disp([' '])
         
      end % End if echo_flag
 
      
      

   case 50 % ---- GLOBAL SENSITIVITY ANALYSIS ------------------------------------------------------

      if echo_flag
         % Clear screen and display message
         disp([' '])
         disp('GLOBAL SENSITIVITY ANALYSIS is running, please wait... (Ctrl+C breaks)')
      end

      % Run the analysis
      [ svrdata, sobolresults, probdata ] = Sobol_SA(1,probdata,analysisopt,gfundata,femodel,randomfield);

      % Store sobolresults in analysisopt structure
      analysisopt.sobolresults = sobolresults;

      if echo_flag
         
         % Display results
         
         if isfield(analysisopt,'SVR')
            
            optSVR     = analysisopt.SVR;
            SVRbasis   = analysisopt.SVRbasis;
            SVR_Nbasis = analysisopt.SVR_Nbasis;
            
            if strcmp(optSVR ,'yes')
               
               disp([' '])
               disp(['..............................................................................................'])
               disp([' '])
               disp(['RESULTS FROM RUNNING SVR-BASED GLOBAL SENSITIVITY ANALYSIS' ])
               disp([' '])
               disp('The following parameters are now available in your current workspace:')
               disp('   sobolresults.First       = First order indices')
               disp('   sobolresults.FirstSecond = First order, second order indices, ...')
               disp('   sobolresults.Total       = Total indices')
               disp('   sobolresults.All         = All indices')
               disp('   sobolresults.nfun        = Number of calls to the limit-state function')
               disp('   svrdata.hyperparameters  = SVR hyperparameters')
               disp('   svrdata.error            = Least square error')
               disp(['..............................................................................................'])

               if analysisopt.first_indices == 1 & analysisopt.all_indices == 0 & analysisopt.total_indices == 0
                  disp('First order indices:'), disp(sobolresults.First)
               end
               if analysisopt.first_indices == 1 & analysisopt.all_indices == 1  & analysisopt.total_indices == 0
                  disp('First order indices:'), disp(sobolresults.FirstSecond(1:size(probdata.marg,1),:))
                  disp('All order indices...:'), disp(sobolresults.FirstSecond(size(probdata.marg,1)+1:size(sobolresults.FirstSecond,1) ,:))
               end
               if analysisopt.first_indices == 1 & analysisopt.all_indices == 0 & analysisopt.total_indices == 1
                  disp('First order indices:'), disp(sobolresults.First)
                  disp('Total indices:'), disp(sobolresults.Total)
                  % disp('All indices:'), disp(sobolresults.All)
               end
               if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
                  disp('First order indices:'), disp(sobolresults.FirstSecond(1:size(probdata.marg,1),:))
                  disp('All order indices...:'), disp(sobolresults.FirstSecond(size(probdata.marg,1)+1:size(sobolresults.FirstSecond,1) ,:))
                  disp('Total indices:'), disp(sobolresults.Total)
                  % disp('All indices:'), disp(sobolresults.All)
               end
               if analysisopt.total_indices == 1 & analysisopt.first_indices == 0 & analysisopt.all_indices == 0
                  disp('Total indices:'), disp(sobolresults.Total)
               end
               disp('Number of calls to the limit-state function:'), disp(sobolresults.nfun)
               disp(['..............................................................................................'])
               disp([' '])

            else % Else if strcmp(optSVR ,'yes')

               disp([' '])
               disp(['..............................................................................................'])
               disp([' '])
               disp(['RESULTS FROM RUNNING GLOBAL SENSITIVITY ANALYSIS' ])
               disp([' '])
               disp('The following parameters are now available in your current workspace:')
               disp('   sobolresults.First       = First order indices')
               disp('   sobolresults.FirstSecond = First order, second order indices, ...')
               disp('   sobolresults.Total       = Total indices')
               disp('   sobolresults.All         = All indices')
               disp('   sobolresults.nfun        = Number of calls to the limit-state function')
               disp(['..............................................................................................'])

               if analysisopt.first_indices == 1 & analysisopt.all_indices == 0 & analysisopt.total_indices == 0
                  disp('First order indices:'), disp(sobolresults.First)
               end
               if analysisopt.first_indices == 1 & analysisopt.all_indices == 1  & analysisopt.total_indices == 0
                  disp('First order indices:'), disp(sobolresults.FirstSecond(1:size(probdata.marg,1),:))
                  disp('All order indices...:'), disp(sobolresults.FirstSecond(size(probdata.marg,1)+1:size(sobolresults.FirstSecond,1) ,:))
               end
               if analysisopt.first_indices == 1 & analysisopt.all_indices == 0 & analysisopt.total_indices == 1
                  disp('First order indices:'), disp(sobolresults.First)
                  disp('Total indices:'), disp(sobolresults.Total)
                  % disp('All indices:'), disp(sobolresults.All)
               end
               if analysisopt.first_indices == 1 & analysisopt.all_indices == 1 & analysisopt.total_indices == 1
                  disp('First order indices:'), disp(sobolresults.FirstSecond(1:size(probdata.marg,1),:))
                  disp('All order indices...:'), disp(sobolresults.FirstSecond(size(probdata.marg,1)+1:size(sobolresults.FirstSecond,1) ,:))
                  disp('Total indices:'), disp(sobolresults.Total)
                  % disp('All indices:'), disp(sobolresults.All)
               end
               if analysisopt.total_indices == 1 & analysisopt.first_indices == 0 & analysisopt.all_indices == 0
                  disp('Total indices:'), disp(sobolresults.Total)
               end
               disp(' Number of calls to the limit-state function:'), disp(sobolresults.nfun)
               disp('( Warning: nfun accounts for the number of replications NbCal )')
               disp(['..............................................................................................'])
               disp([' '])
               
            end % End if strcmp(optSVR ,'yes')
            
         end % End if isfield(analysisopt,'SVR')
         
      end % End if echo_flag




   otherwise % --------------------------------------------------------------------------------------
      
      disp(' ');
      disp('  An invalid choice was entered.');
      disp(' ');

      
   
   
end % --------------------------------------------------------------------------------------


if echo_flag
   toc
end
