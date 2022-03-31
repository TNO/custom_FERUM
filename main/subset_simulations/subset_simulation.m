function [ subsetsimulationresults, probdata ] = subset_simulation(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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


global nfun


if isfield(analysisopt,'ss_restart_from_step')
    ss_restart_from_step = analysisopt.ss_restart_from_step;
else
    ss_restart_from_step = -inf;
end

% Load data from previous subset step, if necessary ( ss_restart_from_step >= 0 )
if ss_restart_from_step >= 0
    eval(['load ss_restart_from_step_' num2str(ss_restart_from_step) '.mat']);
    twister('state',SubsetData.Seed(:,ss_restart_from_step+2)');
end


if ss_restart_from_step < 0

    nfun = 0;

    if ~isfield(analysisopt,'flag_plot')
        analysisopt.flag_plot = 0;
    end
    flag_plot = analysisopt.flag_plot;
    if ~isfield(analysisopt,'flag_plot_gen')
        analysisopt.flag_plot_gen = 0;
    end
    flag_plot_gen = analysisopt.flag_plot_gen;
    if ~isfield(analysisopt,'flag_cov_pf_bounds')
        analysisopt.flag_cov_pf_bounds = 1;
    end
    flag_cov_pf_bounds = analysisopt.flag_cov_pf_bounds;

    % Extract model data
    marg = probdata.marg;
    R = probdata.correlation;
    transf_type = probdata.transf_type;

    % Find number of random variables
    nrv = size(marg,1);

    if isfield(analysisopt,'echo_flag')
        echo_flag = analysisopt.echo_flag;
    else
        echo_flag = 1;
    end

    block_size = analysisopt.block_size;

    rand_generator = analysisopt.rand_generator;
    stdv1 = analysisopt.stdv_sim;
    num_sim = analysisopt.num_sim;

    if isfield(analysisopt,'pf_target')
        pf_target = analysisopt.pf_target;
    else
        analysisopt.pf_target = 0.1;
    end
    SubsetData.pf_target = analysisopt.pf_target;

    if isfield(analysisopt,'width')
        width = analysisopt.width;
    else
        analysisopt.width = 2;
    end
    SubsetData.width = analysisopt.width;

    point = zeros(nrv,1);

    % Modify correlation matrix and perform Cholesky decomposition
    if ~isfield(probdata,'Lo')

        if transf_type == 3

            % Compute corrected correlation coefficients
            switch probdata.Ro_method
                case 0
                    Ro = mod_corr( marg, R );
                case 1
                    [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , 0);
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

    % Establish covariance matrix, its Cholesky decomposition, and its inverse
    covariance = stdv1^2 * eye(nrv);
    chol_covariance = stdv1 * eye(nrv);  % chol_covariance = chol(covariance);
    inv_covariance = 1/stdv1 * eye(nrv); % inv_covariance = inv(covariance);

end


% Subset Simulation - Step 0 (Monte Carlo Simulation)

SubsetData.Nb_step = 0; nb_step = 0;
%asl SubsetData.Seed = twister('state')';
%asl
currgensetting  = rng;                  % Save the current generator settings (the default settings are the Mersenne Twister with seed 0)
SubsetData.Seed = currgensetting.State';
%asl

if ss_restart_from_step == -1
    save ss_restart_from_step_0.mat SubsetData
end

if ss_restart_from_step < 0

    % Initializations
    k = 0;
    percent_done = 0;

    SubsetData.U = zeros(nrv,num_sim);
    SubsetData.G = zeros(1,num_sim);

end

if nrv == 2 & exist('flag_plot') == 1 & flag_plot == 1
    close all
end

if ss_restart_from_step < 0

    while k < num_sim

        block_size = min( block_size, num_sim-k );
        k = k + block_size;

        % Generate realizations of random U-vector
        switch rand_generator
            case 0
                allu = point*ones(1,block_size) + chol_covariance * randn(nrv,block_size);
            otherwise
                % modern MATLAB has Mersenne Twister for generating pseudo
                % random numbers so we use that (FERUM's Twister throws an
                % error), `rand`'s default is the Twister algorithm
                allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(rand(nrv,block_size));
        end

        SubsetData.U(:,(k-block_size+1):k) = allu;

        % Transform into original space
        allx = u_to_x(allu,probdata);

        % Evaluate limit-state function
        [ allg, dummy ] = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
        SubsetData.G((k-block_size+1):k) = allg;

        if floor( k/num_sim * 20 ) > percent_done
            percent_done = floor( k/num_sim * 20 );
            if echo_flag
                fprintf(1,'Subset step #%d - %d%% complete\n',SubsetData.Nb_step,percent_done*5)
            end
        end

    end

    SubsetData.Neval   = num_sim;
    SubsetData.N       = num_sim;
    SubsetData.Indices = [1 num_sim];
    SubsetData.y       = [];
    SubsetData.Indgerm = [];
    SubsetData.AccRate = [];
    SubsetData.p       = [];
    if isfield(analysisopt,'flag_cov_pf_bounds') & analysisopt.flag_cov_pf_bounds == 1
        SubsetData.cov_pf_step  = [];
    end


    % Find y-threshold value.
    % Carry out intermediate calculations for the final estimation
    % of the coefficient variation of the failure probability pf.
    SubsetData = subset_y_threshold(SubsetData,analysisopt);

    %asl SubsetData.Seed = [ SubsetData.Seed twister('state')' ];
    %asl
    currgensetting  = rng;                  % Save the current generator settings (the default settings are the Mersenne Twister with seed 0)
    SubsetData.Seed = [ SubsetData.Seed currgensetting.State' ];
    %asl

    if ss_restart_from_step == -1
        ss_restart_from_step0 = ss_restart_from_step;
        analysisopt0 = analysisopt;
        clear ss_restart_from_step
        analysisopt = rmfield(analysisopt,'ss_restart_from_step');
        save ss_restart_from_step_0.mat
        ss_restart_from_step = ss_restart_from_step0;
        analysisopt = analysisopt0;
    end

end


if nrv == 2 & exist('flag_plot') == 1 & flag_plot == 1
    SubsetData = subset_subset_graph(SubsetData,0);
end


if SubsetData.y ~= 0

    % Subset Simulation - Step 1

    SubsetData.Nb_step = SubsetData.Nb_step + 1; nb_step = SubsetData.Nb_step;

    if SubsetData.Nb_step > ss_restart_from_step
        if ~( ss_restart_from_step < 0 )
            analysisopt.ss_restart_from_step = -1;
            ss_restart_from_step = -1;
        end
    end

    if ss_restart_from_step < 0
        % Subset Simulation - Step 1 - Run simulations
        SubsetData = subset_subset_step(SubsetData,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
    end

    if nrv == 2 & exist('flag_plot') == 1 & flag_plot == 1
        SubsetData = subset_subset_graph(SubsetData,1);
    end

end

% Loop until y-threshold equal to zero
while SubsetData.y(end) ~= 0

    if ss_restart_from_step < 0
        % Find y-threshold value.
        % Carry out intermediate calculations for the final estimation
        % of the coefficient variation of the failure probability pf.
        SubsetData = subset_y_threshold(SubsetData,analysisopt);
    end

    %asl SubsetData.Seed = [ SubsetData.Seed twister('state')' ];
    %asl
    currgensetting  = rng;                  % Save the current generator settings (the default settings are the Mersenne Twister with seed 0)
    SubsetData.Seed = [ SubsetData.Seed currgensetting.State' ];
    %asl

    if ss_restart_from_step == -1
        ss_restart_from_step0 = ss_restart_from_step;
        analysisopt0 = analysisopt;
        clear ss_restart_from_step
        analysisopt = rmfield(analysisopt,'ss_restart_from_step');
        eval(['save ss_restart_from_step_' num2str(SubsetData.Nb_step) '.mat']);
        ss_restart_from_step = ss_restart_from_step0;
        analysisopt = analysisopt0;
    end

    if nrv == 2 & exist('flag_plot') == 1 & flag_plot == 1
        SubsetData = subset_subset_graph(SubsetData,0);
    end

    if SubsetData.y(end) == 0, break, end

    % Subset Simulation - Step > 1
    SubsetData.Nb_step = SubsetData.Nb_step + 1; nb_step = SubsetData.Nb_step;

    if SubsetData.Nb_step > ss_restart_from_step
        if ~( ss_restart_from_step < 0 )
            analysisopt.ss_restart_from_step = -1;
            ss_restart_from_step = -1;
        end
    end

    if ss_restart_from_step < 0
        % Subset Simulation - Step > 1 - Run simulations
        SubsetData = subset_subset_step(SubsetData,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
    end

    if nrv == 2 & exist('flag_plot') == 1 & flag_plot == 1 & SubsetData.y(end) ~= 0
        SubsetData = subset_subset_graph(SubsetData,1);
    end

end

% Assess approximation bounds for the coeficient of variation of the failure probability pf
if isfield(analysisopt,'flag_cov_pf_bounds') & analysisopt.flag_cov_pf_bounds == 1
    SubsetData.cov_pf_bounds = [ sqrt(sum((SubsetData.cov_pf_step).^2)) sqrt(sum(sum((SubsetData.cov_pf_step)'*(SubsetData.cov_pf_step)))) ];
end

if nrv == 2 & exist('flag_plot') == 1 & flag_plot == 1
    SubsetData = subset_subset_graph(SubsetData,3);
end

% Plot generalized reliability index vs. threshold value, with +/-2 stdv interval, for each subset step.
% Normal and lognormal hypothesis based on first subset step (CMC) samples
if echo_flag
    cov_pf_bounds = [];
    Pf = [];
    Beta = []; Beta_inf = []; Beta_sup = [];
    for nb_step = 0:SubsetData.Nb_step
        cov_pf_inf = sqrt(sum((SubsetData.cov_pf_step(1:(nb_step+1))).^2));
        cov_pf_sup = sqrt(sum(sum((SubsetData.cov_pf_step(1:(nb_step+1)))'*(SubsetData.cov_pf_step(1:(nb_step+1))))));
        cov_pf_bounds = [ cov_pf_bounds [ cov_pf_inf ; cov_pf_sup ] ];
        pf = prod(SubsetData.p(1:(nb_step+1)));
        Pf = [ Pf pf ];
        Beta  = [ Beta -inv_norm_cdf(pf) ];
        Beta_inf = [ Beta_inf -inv_norm_cdf(pf+2*cov_pf_inf*pf) ];
        Beta_sup = [ Beta_sup -inv_norm_cdf(pf-2*cov_pf_inf*pf) ];
        if nb_step == 0
            mean0 = mean(SubsetData.G(1:SubsetData.Indices(1,2)));
            stdv0 = std(SubsetData.G(1:SubsetData.Indices(1,2)));
          	cov0 = stdv0/mean0;
          	zeta0 = (log(1+cov0^2))^0.5;
          	lambda0 = log(mean0) - 0.5*zeta0^2;
        end
    end
    [ Pf-2*cov_pf_bounds(1,:).*Pf ; Pf; Pf+2*cov_pf_bounds(1,:).*Pf ];
    [ Beta_sup; Beta; Beta_inf ];
    figure
    hbeta = errorbar(SubsetData.y,Beta,Beta-Beta_inf,Beta_sup-Beta);
    hold on
    hnorm = plot(SubsetData.y,-(SubsetData.y-mean0)/stdv0,'r--');
    hlogn = plot(SubsetData.y,-(log(SubsetData.y)-lambda0)/zeta0,'r:');
    set(hbeta,'LineWidth',2)
    set(hnorm,'LineWidth',2)
    set(hlogn,'LineWidth',2)
    set(gca,'FontSize',14);
    xlabel('{\ity_i}','FontSize',14);
    ylabel('{\it\beta_i}','FontSize',14);
end


SubsetData.pf = prod(SubsetData.p);


subsetsimulationresults.pf         = SubsetData.pf;
subsetsimulationresults.beta       = -inv_norm_cdf(SubsetData.pf);
subsetsimulationresults.nfun       = nfun;
subsetsimulationresults.SubsetData = SubsetData;

if isfield(analysisopt,'flag_cov_pf_bounds') & analysisopt.flag_cov_pf_bounds == 1
    subsetsimulationresults.cov_pf  = SubsetData.cov_pf_bounds(1);
end
