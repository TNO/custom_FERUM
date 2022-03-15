function P = ferum_cdf(type,x,param)

% Cumulative Density Function
%
%!   P = ferum_cdf(type,x,param)
%
%   Evaluates the cumulative distribution function and returns the probability.
%
%   Output: - P           = cumulative density value
%   Input:  - type        = probability distribution type (1: normal, 2: lognormal, ...)
%           - x           = 'Abscissa' value(s)
%!          - param(1) = parameter #1 of the random variable
%!          - param(2) = parameter #2 of the random variable (if applicable)
%!          - param(3) = parameter #3 of the random variable (if applicable)
%!          - param(4) = parameter #4 of the random variable (if applicable)

switch type
    
    case 1  % Normal marginal distribution
        
        mean    = param(1);
        stdv    = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = (0.5+erf(((x-mean)/stdv)/sqrt(2))/2).^n;
        
    case 2  % Lognormal marginal distribution
        
        lambda  = param(1);
        zeta    = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        z       = ( log(x) - lambda ) / zeta;
        P       = (0.5+erf(z/sqrt(2))/2).^n;
        
    case 3  % Gamma distribution
        
        lambda  = param(1);
        k       = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = (gammainc(lambda*x,k)).^n;
        
    case 4  % Shifted exponential distribution
        
        lambda  = param(1);
        x_zero  = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = (1 - exp( -lambda*(x-x_zero) )).^n;
        
    case 5  % Shifted Rayleigh distribution
        
        a       = param(1);
        x_zero  = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = (1 - exp( -0.5*( (x-x_zero)/a ).^2 )).^n;
        
    case 6  % Uniform distribution
        
        a       = param(1);
        b       = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = ((x-a) / (b-a)).^n;
        
    case 7  % Beta distribution
        
        q       = param(1);
        r       = param(2);
        a       = param(3);
        b       = param(4);
        
        x01     = (x-a) / (b-a);
        P       = betainc(x01,q,r);
        
    case 8  % Chi-square distribution
        
        nu      = param(1);
        lambda  = 0.5;
        k       = nu/2;
        %!param(2) not relevant
        %!param(3) not relevant
        n       = param(4);
        
        P       = (gammainc(lambda*x,k)).^n;
        
    case 11  % Type I largest value distribution ( same as Gumbel distribution )
        
        u_n     = param(1);    % location
        a_n     = param(2);    % inverse scale
        %!param(3) not relevant
        n       = param(4);
        
        P       = (exp( -exp( -a_n*(x-u_n) ) )).^n;
        
    case 12  % Type I smallest value distribution
        
        u_1     = param(1);
        a_1     = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = (1 - exp( -exp( a_1*(x-u_1) ) )).^n;
        
    case 13  % Type II largest value distribution
        
        u_n     = param(1);
        k       = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P = (exp( -(u_n./x).^k )).^n;
        
    case 14  % Type III smallest value distribution
        
        u_1     = param(1);
        k       = param(2);
        epsilon = param(3);
        n       = param(4);
        
        P       = (1 - exp( -( (x-epsilon) / (u_1-epsilon) ).^k )).^n;
        
    case 15  % Gumbel distribution ( same as type I largest value distribution )
        
        u_n     = param(1);
        a_n     = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = (exp( -exp( -a_n*(x-u_n) ) )).^n;
        
    case 16  % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
        
        u_1     = param(1);
        k       = param(2);
        %!param(3) not relevant
        n       = param(4);
        
        P       = (1 - exp( -(x/u_1).^k )).^n;
        
    case 18  % (Reserved for Laplace marginal distribution)
        
    case 19  % (Reserved for Pareto marginal distribution)
        
    case 20 % Generalized extreme value (GEV) distribution
        k       = param(1); % shape
        sigma   = param(2); % scale
        mu      = param(3); % location
        n       = param(4);
        
        P       = (gevcdf(x, k, sigma, mu)).^n;
        
    case 25 % Three-parameter lognormal (LN3) distribution
        shape   = param(1);
        scale   = param(2);
        thres   = param(3);
        n       = param(4);
        
        P       = (lognorm3cdf(x, shape, scale, thres)).^n;    
        
    case 30 % sample based custom distibution - typically for non-parametric distributions
        % powered option is valid, its effect is directly (numerically) incorporated into the kernel distribution (kernel_sample_*.mat)
        ID      = param(1);
        %!param(2) not relevant
        %!param(3) not relevant
        %!param(4) not relevant
        
        P       = custom_cdf(x, ID, 'sample');
        
    case 31 % vector based custom distibution - typically for non-parametric distributions
        % USE 32 INSTEAD OF THIS!        
        % powered option is not valid
        ID      = param(1);
        %!param(2) not relevant
        %!param(3) not relevant
        %!param(4) not relevant
        
        P       = custom_cdf(x, ID, 'point');
    
    case 32 % generalized non-parametric distribution - vector based
        % USE THIS INSTEAD OF 31!
        % powered option is not valid
        ID          = param(1);
        shift       = param(2);
        scale       = param(3);
        %!param(4) not relevant
        
        P           = nonparametric_cdf(x, ID, shift, scale);
    case 33 % hardcoded
        % USE THIS INSTEAD OF 31!
        % powered option is not valid
        ID          = param(1);
        shift       = param(2);
        scale       = param(3);
        %!param(4) not relevant
        
        P           = hardcoded_cdf(x, ID, shift, scale);

    case 51  % Truncated normal marginal distribution
        
        mean    = param(1);
        stdv    = param(2);
        xmin    = param(3);
        xmax    = param(4);
        
        P       = ( erf(((x-mean)/stdv)/sqrt(2))/2 - erf(((xmin-mean)/stdv)/sqrt(2))/2 ) / ...
                  ( erf(((xmax-mean)/stdv)/sqrt(2))/2 - erf(((xmin-mean)/stdv)/sqrt(2))/2 );
        Imin    = find(x<xmin); P(Imin) = 0;
        Imax    = find(x>xmax); P(Imax) = 1;
        
    otherwise
        
end
