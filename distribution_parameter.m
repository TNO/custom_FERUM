function marg  = distribution_parameter(marg)

% This function computes:
%    - the parameters corresponding to the marginal distributions knowing the mean and stdv
%                                            or
%    - the mean and stdv corresponding to the marginal distributions knowing the parameters
%
% CALLED CUSTOM FUNCTION(S):
% kde. m    http://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation

% clear potential persistent variables
clear custom_pdf custom_cdf custom_invcdf
clear nonparametric_pdf nonparametric_cdf nonparametric_invcdf
clear lognlognormcdf

nrv = size(marg,1);
flag_31 = 0;
flag_32 = 0;

for i=1:nrv
    
    switch marg(i,1)
        
        case 1 % Normal distribution
            
            if marg(i,9) == 1
                mean = marg(i,5);
                stdv = marg(i,6);
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                marg(i,5) = mean;
                marg(i,6) = stdv;
                marg(i,7) = 0;
                %!marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 2 % Lognormal distribution
            
            if marg(i,9) == 1
                lambda = marg(i,5);
                zeta = marg(i,6);
                mean = exp(lambda+0.5*(zeta^2));
                stdv = exp(lambda+0.5*(zeta^2)) * (exp(zeta^2)-1)^0.5;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                cov = stdv/mean;
                zeta = (log(1+cov^2))^0.5;
                lambda = log(mean) - 0.5*zeta^2;
                marg(i,5) = lambda;
                marg(i,6) = zeta;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 3 % Gamma distribution
            
            if marg(i,9) == 1
                lambda = marg(i,5);
                k      = marg(i,6);
                mean = k/lambda;
                stdv = (k^0.5)/lambda;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                lambda = mean / stdv^2 ;
                k = mean^2 / stdv^2;
                marg(i,5) = lambda;
                marg(i,6) = k;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 4 % Shifted exponential distribution
            
            if marg(i,9) == 1
                lambda = marg(i,5);
                x_zero = marg(i,6);
                mean = x_zero + 1/lambda;
                stdv = 1/lambda;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                x_zero= mean - stdv;
                lambda= 1/stdv;
                marg(i,5) = lambda;
                marg(i,6) = x_zero;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 5 % Shifted Rayleigh distribution
            
            if marg(i,9) == 1
                a      = marg(i,5);
                x_zero = marg(i,6);
                mean = x_zero + a*(pi/2)^0.5;
                stdv = a*(2-pi/2)^0.5;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                a = stdv / (2-pi/2)^0.5 ;
                x_zero = mean - (pi/(4-pi))^0.5 * stdv;
                marg(i,5) = a;
                marg(i,6) = x_zero;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 6 % Uniform distribution
            
            if marg(i,9) == 1
                a = marg(i,5);
                b = marg(i,6);
                mean = (a+b)/2 ;
                stdv = (b-a)/(2*(3)^0.5);
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                a = mean - 3^0.5 * stdv;
                b = mean + 3^0.5 * stdv;
                marg(i,5) = a;
                marg(i,6) = b;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 7 % Beta distribution
            
            if marg(i,9) == 1
                q = marg(i,5);
                r = marg(i,6);
                a = marg(i,7);
                b = marg(i,8);
                mean = a + q*(b-a)/(q+r);
                stdv = ((b-a)/(q+r))*(q*r/(q+r+1))^0.5;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                a    = marg(i,7);
                b    = marg(i,8);
                par_guess = 1;
                par = fminsearch('betapar',par_guess,optimset('fminsearch'),a,b,mean,stdv);
                q = par;
                r = q*(b-a)/(mean-a) - q;
                marg(i,5) = q;
                marg(i,6) = r;
                marg(i,7) = a;
                marg(i,8) = b;
            end
            
        case 8 % Chi-square distribution
            
            if marg(i,9) == 1
                lambda = 0.5;
                nu     = marg(i,5);
                mean = nu/(2*lambda);
                stdv = ((nu/2)^0.5)/lambda;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                lambda = 0.5;
                mean_test = lambda*stdv^2;
                if mean/mean_test < 0.95 || mean/mean_test > 1.05
                    error('Error when using Chi-square distribution. Mean and stdv should be given such that mean = 0.5*stdv.^2\n')
                end
                nu = 2*(mean^2/stdv^2);
                marg(i,5) = nu;
                marg(i,6) = 0;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 11 % Type I largest value distribution ( same as Gumbel distribution )
            
            am = 0.57721566490153286060651209008240243104215933593992;
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
            if marg(i,9) == 1
                u_n = marg(i,5);
                a_n = marg(i,6);
                mean = u_n + am/a_n;
                stdv = pi/(a_n*6^0.5);
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                n       = marg(i,8);
                mean    = marg(i,2);
                stdv    = marg(i,3);
                a_n     = pi/(stdv*sqrt(6));
                u_n     = mean - (am*stdv*sqrt(6))/pi;
                
                % WARNING!
                % set the starting point for FORM to updated mean
                if isnan(marg(i,4))
                    u_n_n   = u_n + 1/a_n*log(n);
                    a_n_n   = a_n;
                    mean_n  = u_n_n + am/a_n_n;                
                    marg(i,4) = mean_n;
                end
                marg(i,5) = u_n;
                marg(i,6) = a_n;
                marg(i,7) = 0;
            end
            

            
        case 12 % Type I smallest value distribution
            
            if marg(i,9) == 1
                u_1 = marg(i,5);
                a_1 = marg(i,6);
                mean = u_1 - 0.5772156649/a_1 ;
                stdv = pi/(a_1*6^0.5);
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                a_1 = pi/(stdv*sqrt(6));
                u_1 = mean + (0.5772156649*stdv*sqrt(6))/pi;
                marg(i,5) = u_1;
                marg(i,6) = a_1;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 13 % Type II largest value distribution
            
            if marg(i,9) == 1
                u_n = marg(i,5);
                k   = marg(i,6);
                mean = u_n*gamma(1-1/k);
                stdv = u_n*(gamma(1-2/k)-(gamma(1-1/k))^2)^0.5;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                par_guess = [2.000001 10e+07];
                par = fzero('typIIlargestpar',par_guess,optimset('fzero'),mean,stdv);
                k = par;
                u_n = mean/gamma(1-1/k);
                marg(i,5) = u_n;
                marg(i,6) = k;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 14 % Type III smallest value distribution
            
            if marg(i,9) == 1
                u_1     = marg(i,5);
                k       = marg(i,6);
                epsilon = marg(i,7);
                mean = epsilon + (u_1-epsilon)*gamma(1+1/k);
                stdv = (u_1-epsilon)*(gamma(1+2/k)-gamma(1+1/k)^2)^0.5;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean    = marg(i,2);
                stdv    = marg(i,3);
                epsilon = marg(i,7);
                meaneps = mean - epsilon;
                par_guess = [0.1 10e+07];
                par = fzero('typIIIsmallestpar',par_guess,optimset('fzero'),meaneps,stdv);
                k = par;
                u_1 = meaneps/gamma(1+1/k)+epsilon;
                marg(i,5) = u_1;
                marg(i,6) = k;
                marg(i,7) = epsilon;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 15 % Gumbel distribution ( same as type I largest value distribution )
            
            if marg(i,9) == 1
                u_n = marg(i,5);
                a_n = marg(i,6);
                mean = u_n + 0.5772156649/a_n ;
                stdv = pi/(a_n*6^0.5);
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9)== 0
                mean = marg(i,2);
                stdv = marg(i,3);
                a_n = pi/(stdv*sqrt(6));
                u_n = mean - (0.5772156649*stdv*sqrt(6))/pi;
                marg(i,5) = u_n;
                marg(i,6) = a_n;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
            
            if marg(i,9) == 1
                u_1 = marg(i,5);
                k = marg(i,6);
                mean = u_1*gamma(1+1/k);
                stdv = u_1*(gamma(1+2/k)-gamma(1+1/k)^2)^0.5;
                marg(i,2) = mean;
                marg(i,3) = stdv;
            elseif marg(i,9) == 0
                mean = marg(i,2);
                stdv = marg(i,3);
                epsilon = 0;
                meaneps = mean - epsilon;
                par_guess = [0.1 10e+07];
                par = fzero('typIIIsmallestpar',par_guess,optimset('fzero'),meaneps,stdv);
                k = par;
                u_1 = meaneps/gamma(1+1/k)+epsilon;
                marg(i,5) = u_1;
                marg(i,6) = k;
                marg(i,7) = 0;
                %marg(i,8) = 0;
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 18 % (Reserved for Laplace distribution)
            
        case 19 % (Reserved for Pareto distribution)
            
        case 20 % Generalized extreme value (GEV) distribution
            if marg(i,9) == 1
                k       = marg(i,5);
                sigma   = marg(i,6);
                mu      = marg(i,7);
                [mean, var] = gevstat(k, sigma, mu);
                marg(i,2) = mean;
                marg(i,3) = sqrt(var);
                if isnan(marg(i,4))
                    marg(i,4) = mean;
                end
            elseif marg(i,9) == 0
                error('Generalized extreme value distribution (3 parameters) can be defined only with 3 distribution parameters! Definition with mean and standard deviation is ambiguous.')
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 25 % Three-parameter lognormal (LN3) distribution
            if marg(i,9) == 1
                shape       = marg(i,5);
                scale       = marg(i,6);
                thres       = marg(i,7);
                mom         = lognorm3stat(shape, scale, thres, 'par');
                marg(i,2)   = mom(1);
                marg(i,3)   = mom(1)*mom(2);
                if isnan(marg(i,4))
                    marg(i,4) = mean;
                end
            elseif marg(i,9) == 0
                error('Three-parameter lognormal distribution (3 parameters) can be defined only with 3 distribution parameters! Definition with mean and standard deviation is ambiguous.')
            end
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            
        case 30 % Custom distibution, defined by a vector of sample
            % Construct the kernel pdf and cdf functions and save it to temp/kernel_dist_*.mat
            
            %! set the powered distribution option
            if isnan(marg(i,8)) || marg(i,8) <= 0
                marg(i,8) = 1;
            end
            % value of marg(i,9) does not matter
            ID       = marg(i,5);
            n        = marg(i,8);
            
            load(['tmp\sample_',num2str(ID),'.mat'], 'sample')
            %kernel density implementation from file exchange (Botev)
            % with default 2^14 meshpoints and min  max (range) +0.1*range
            [~,pdf,x_grid,cdf] = kde(sample);
            
            x_grid = x_grid';
            % apply the power, makes difference only if it is provided and differs from zero or NaN
            cdf = cdf.^n;
            %...............................................................
            % Correct the kernel cdf to make strictly monotonic increasing, this should affect only the regions not relevant to the analysis
            % I am not satisfied with this solution
            th_l = 10^-10;
            th_r = 10^-10;
            
            while any(diff(cdf)<=0);
                
                % correct the left tail
                %any(diff(cdf((cdf < th_l)))<=0)
                leftflag = any(diff(cdf((cdf < 0.5)))<=0);
                if leftflag == 1
                    idx = sum(cdf < th_l);
                    cdf(cdf < th_l) = linspace(cdf(1), cdf(idx), idx);
                    th_l = 5*th_l;
                end
                
                % correct the right tail
                %any(diff(cdf((cdf >= 1-th_r)))<=0)
                rightflag = any(diff(cdf((cdf >= 0.5)))<=0);
                if rightflag == 1
                    idx = length(cdf) - sum(cdf >= 1-th_r) + 1;
                    cdf(cdf >= 1-th_r) = linspace(cdf(idx), cdf(end), sum(cdf >= 1-th_r));
                    th_r = 5*th_r;
                end
            end
            %            disp(['Kernel correction left threshold:', num2str(th_l)])
            %            disp(['Kernel correction right threshold:', num2str(th_r)])
            %...............................................................
            
            % properties of the kernel density estimation, should be compared to the design point to see how 'reliable' is the solution
            prop.samplemax = max(sample);
            prop.samplemin = min(sample);
            prop.thresholdlower = th_l;
            prop.thresholdupper = 1-th_r;
            
            if n == 1 % not powered distribution, to avoid finite difference calculation
                % do nothing, use the pdf provided by kde()
            else % powered distribution
                Fn  = @(x) interp1(x_grid, cdf, x);
                
                pdf = cfd(Fn, x_grid); % might be moved to ferum_pdf to reduce computation demand (not significant)
                
                % correct the ends if necessary
                if isnan(pdf(1))
                    pdf(1) = pdf(2)-1e-10;
                end
                if isnan(pdf(end))
                    pdf(end) = pdf(end-1)+1e-10;
                end
                if isnan(pdf)
                    error('There is some problem with the numerical derivation of the powered cdf function!')
                end
            end
            
            % set the mean and std based on the sample and power
            dxc         = diffc(x_grid);
            mean_i      = sum(pdf.*x_grid.*dxc);
            stdv_i      = sqrt(sum(pdf.*(x_grid-mean_i).^2.*dxc));
            marg(i,2)   = mean_i;
            marg(i,3)   = stdv_i;
            
            % set the starting point to mean if not specified
            if isnan(marg(i,4))
                marg(i,3) = mean_i;
            end
            
            
            save(['tmp\kernel_sample_', num2str(ID), '.mat'], 'x_grid', 'pdf', 'cdf', 'prop')
            
        case 31 % Custom distibution, defined by vector of points!
            % USE 32 INSTEAD OF THIS!
            % powered option is not valid (not implemented)
            ID       = marg(i,5);
            
            load(['tmp\vector_distr_', num2str(ID), '.mat'], 'x_grid', 'pdf', 'cdf')
            
            % =================================================
            % tmp solution!
            idx         = isnan(x_grid) | isnan(pdf) | isnan(cdf);
            x_grid      = x_grid(~idx);
            pdf         = pdf(~idx);
            cdf         = cdf(~idx);
            % =================================================
            
            mean_i      = trapz(x_grid, pdf.*x_grid);
            stdv_i      = sqrt(trapz(x_grid, pdf.*(x_grid-mean_i).^2));
            marg(i,2)   = mean_i;
            marg(i,3)   = stdv_i;
            
            % set the starting point to mean if not specified
            if isnan(marg(i,4))
                marg(i,4) = mean_i;
            end
            
            flag1 = interp1(cdf, x_grid, normcdf(-6));
            flag2 = interp1(cdf, x_grid, normcdf(6));
            if any(isnan([flag1, flag2]))
               warning('If correlated it won''t work!!') 
            end
            
            S(ID).x_grid   = x_grid;
            S(ID).pdf      = pdf;
            S(ID).cdf      = cdf;
            flag_31        = 1;

         case 32
            % USE THIS INSTEAD OF 31!
            % powered option is not valid (not implemented)
            % shift
            if isnan(marg(i,6))
                marg(i,6) = 0;
            end
            % scale
            if isnan(marg(i,7))
                marg(i,7) = 1;
            end
            
            
            ID          = marg(i,5);
            shift       = marg(i,6);
            scale       = marg(i,7);
            
            load(['tmp\vector_distr_', num2str(ID), '.mat'], 'x_grid', 'pdf', 'cdf')
            
            % =================================================
            % tmp solution!
            idx         = isnan(x_grid) | isnan(pdf) | isnan(cdf);
            x_grid      = x_grid(~idx);
            pdf         = pdf(~idx);
            cdf         = cdf(~idx);
            % =================================================
            
            mean_i      = trapz(x_grid, pdf.*x_grid);
            stdv_i      = sqrt(trapz(x_grid, pdf.*(x_grid-mean_i).^2));
            marg(i,2)   = scale*mean_i + shift;
            marg(i,3)   = scale*stdv_i;
            
            % set the starting point to mean if not specified
            if isnan(marg(i,4))
                marg(i,4) = marg(i,2);
            end
            
            flag1 = interp1(cdf, x_grid, normcdf(-6));
            flag2 = interp1(cdf, x_grid, normcdf(6));
%             if any(isnan([flag1, flag2]))
%                warning('If correlated it won''t work!!') 
%             end
            
            S(ID).x_grid   = x_grid;
            S(ID).pdf      = pdf;
            S(ID).cdf      = cdf;
            flag_32        = 1;
            
        case 33
            if isnan(marg(i,6))
                marg(i,6) = 0;
            end
            % scale
            if isnan(marg(i,7))
                marg(i,7) = 1;
            end
            
            
            ID          = marg(i,5);
            shift       = marg(i,6);
            scale       = marg(i,7);
            

            [~, x_grid, pdf]    = hardcoded_pdf(1, ID);
            [~, ~, cdf]         = hardcoded_cdf(1, ID);
  
            
            mean_i      = trapz(x_grid, pdf.*x_grid);
            stdv_i      = sqrt(trapz(x_grid, pdf.*(x_grid-mean_i).^2));
            marg(i,2)   = scale*mean_i + shift;
            marg(i,3)   = scale*stdv_i;
            
            % set the starting point to mean if not specified
            if isnan(marg(i,4))
                marg(i,4) = marg(i,2);
            end
        case 51 % Truncated normal marginal distribution
            
            if marg(i,9) == 1
                mean = marg(i,5);
                stdv = marg(i,6);
                xmin = marg(i,7);
                xmax = marg(i,8);
                marg(i,2) = mean_norm_truncated(mean,stdv,xmin,xmax);
                marg(i,3) = stdv_norm_truncated(mean,stdv,xmin,xmax);
            end
    end
    
end


for i=1:nrv
    if isnan(marg(i,4))
        marg(i,4) = marg(i,2);
    end
end

if flag_31 || flag_32
   save('tmp\vector_distr_struct.mat', 'S') 
end