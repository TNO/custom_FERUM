function x = ferum_invcdf(type,P,param)

%! Inverse Cumulative Density Function
%
%!   P = ferum_cdf(type,x,param)
%
%   Evaluates the cumulative distribution function and returns the probability.
%
%   Output: - x           = 'Abscissa' value(s)
%   Input:  - type        = probability distribution type (1: normal, 2: lognormal, ...)
%           - P           = cumulative density value
%!          - param(1) = parameter #1 of the random variable
%!          - param(2) = parameter #2 of the random variable (if applicable)
%!          - param(3) = parameter #3 of the random variable (if applicable)
%!          - param(4) = parameter #4 of the random variable (if applicable)


%! the x<0 is not treated! if the problem is well defined the 0 part of cdf never should be reached in the FORM
pon     = param(7);
if pon ~= 1 && any(P < (1-pon))
    error('For an intermittent distribution the search algorithm reached the poff=1-pon region, this indicates error in the probabilistic and/or determinsitic models.')
end

switch type
    
    case 1  % Normal marginal distribution
        
        mean    = param(1);
        stdv    = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = norminv(((P - (1-pon))/pon).^(1/n), mean, stdv);
        %x(i,:) = z(i,:) * marg(i,3) + marg(i,2);
        
    case 2  % Lognormal marginal distribution
        
        lambda  = param(1);
        zeta    = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        mu      = lambda;
        sigma   = zeta;
        
        x       = logninv(((P - (1-pon))/pon).^(1/n), mu, sigma);
        
    case 3  % Gamma distribution
        error('Gamma distribution is not yet fully implemented!')
        %         lambda  = marg(i,5);
        %         k = marg(i,6);
        %         % marg(i,7) not relevant
        %         n       = marg(i,8);
        %         mean = marg(i,2);
        %         for j = 1 : nx
        %             normal_val = normcdf(z(i,j)).^(1/n);
        %             %x(i,j) = fzero('zero_gamma',mean,optimset('fzero'),k,lambda,normal_val); % Doesn't work
        %             x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
        %             %A = k;
        %             %B = -1/(log(lambda^k)*lambda);
        %             %x(i,:) =  gaminv(normcdf(z(i,:)).^(1/n), A, B);
        %             % should be verified
        %         end
        
    case 4  % Shifted exponential distribution
        
        lambda  = param(1);
        x_zero  = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = x_zero + 1/lambda * log( 1 ./ ( 1 - ((P - (1-pon))/pon).^(1/n) ) );
        
    case 5  % Shifted Rayleigh distribution
        
        a       = param(1);
        x_zero  = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = x_zero + a * ( 2*log( 1 ./ (1 - ((P - (1-pon))/pon).^(1/n)) ) ) .^0.5;
        
    case 6  % Uniform distribution
        
        a       = param(1);
        b       = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = a + (b-a) * ((P - (1-pon))/pon).^(1/n);
        
    case 7  % Beta distribution
        error('Beta distribution is not yet fully implemented!')
        %         q       = marg(i,5);
        %         r       = marg(i,6);
        %         a       = marg(i,7);
        %         b       = marg(i,8);
        %         mean    = marg(i,2);
        %         for j = 1 : nx
        %             normal_val = normcdf(z(i,:));
        %             %x01 = fzero('zero_beta',(mean-a)/(b-a),optimset('fzero'),q,r,normal_val); % Doesn't work
        %             x01 = fminbnd('zero_beta',0,1,optimset('fminbnd'),q,r,normal_val);
        %             % Transform x01 from [0,1] to [a,b] interval
        %             x(i,j) = a + x01 * ( b - a );
        %         end
        
    case 8  % Chi-square distribution
        error('Beta distribution is not yet fully implemented!')
%         lambda  = 0.5;
%         nu      = marg(i,5);
%         % marg(i,6) not relevant
%         % marg(i,7) not relevant
%         n       = marg(i,8);
%         k       = nu/2 ;
%         mean    = marg(i,2);
%         for j = 1 : nx
%             normal_val = normcdf(z(i,j)).^(1/n);
%             %x(i,j) = fzero('zero_gamma',mean,optimset('fzero'),k,lambda,normal_val); % Doesn't work
%             x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
%         end
        
    case 11  % Type I largest value distribution ( same as Gumbel distribution )
        
        u_n     = param(1);    % location
        a_n     = param(2);    % inverse scale
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
                                    
        x       = u_n - (1/a_n) * log( log( 1 ./ ((P - (1-pon))/pon).^(1/n)) ) ;
        
    case 12  % Type I smallest value distribution
        
        u_1     = param(1);
        a_1     = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = u_1 + (1/a_1) * log( log( 1 ./ ( 1 - ((P - (1-pon))/pon).^(1/n)) ) ) ;
        
    case 13  % Type II largest value distribution
        
        u_n     = param(1);
        k       = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = u_n * log( 1 ./ ((P - (1-pon))/pon).^(1/n)).^ (-1/k);
        
    case 14  % Type III smallest value distribution
        
        u_1     = param(1);
        k       = param(2);
        epsilon = param(3);
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = epsilon + ( u_1 - epsilon ) * log( 1 ./ ( 1 - ((P - (1-pon))/pon).^(1/n)) ).^(1/k);
        
    case 15  % Gumbel distribution ( same as type I largest value distribution )
        
        u_n     = param(1);
        a_n     = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);

        x       = u_n - (1/a_n) * log( log( 1 ./ ((P - (1-pon))/pon).^(1/n)) );
        
    case 16  % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
        
        u_1     = param(1);
        k       = param(2);
        %!param(3) not relevant
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = u_1 * log( 1 ./ ( 1 - ((P - (1-pon))/pon).^(1/n)) ).^(1/k);
        
    case 18  % (Reserved for Laplace marginal distribution)
        
    case 19  % (Reserved for Pareto marginal distribution)
        
    case 20 % Generalized extreme value (GEV) distribution
        k       = param(1);
        sigma   = param(2);
        mu      = param(3);
        %!param(4) not relevant
        %!param(5) distribution defintion type
        n       = param(6);
        pon     = param(7);
        
        x       = gevinv( ((P - (1-pon))/pon).^(1/n), k, sigma, mu ) ;
        
    case 30 % sample based custom distibution - typically for non-parametric distributions
        % powered option is valid, its effect is directly (numerically) incorporated into the kernel distribution (kernel_sample_*.mat)
        ID      = param(1);
        %!param(2) not relevant
        %!param(3) not relevant
        %!param(4) not relevant
        
% % %         x       = custom_invcdf( P, ID, 'sample');
%         
    error('NYAAAK')
    case 51  % Truncated normal marginal distribution
         error('Beta distribution is not yet fully implemented!')
%         mean    = marg(i,5);
%         stdv    = marg(i,6);
%         xmin    = marg(i,7);
%         xmax    = marg(i,8);
%         x(i,:)  = mean + stdv * inv_norm_cdf( ...
%             normcdf((xmin-mean)/stdv) + ...
%             (normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) * normcdf(z(i,:)) ...
%             );
%         Imin    = find(x(i,:)<xmin); x(i,Imin) = xmin;
%         Imax    = find(x(i,:)>xmax); x(i,Imax) = xmax;
        
    otherwise
        
end
