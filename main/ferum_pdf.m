function p = ferum_pdf(type,x,param)

% Probability Density Function
%
%!   p = ferum_pdf(type,x,param)
%
%   Evaluates probability density function and return the corresponding density value.
%
%   Output: - p           = probability density value
%   Input:  - type        = probability distribution type (1: normal, 2: lognormal, ...)
%           - x           = 'Abscissa' value(s)
%!          - param(1) = parameter #1 of the random variable
%!          - param(2) = parameter #2 of the random variable (if applicable)
%!          - param(3) = parameter #3 of the random variable (if applicable)
%!          - param(4) = parameter #4 of the random variable (if applicable)

switch type
    
    case 1  % Normal distribution
        
        mean = param(1);
        stdv = param(2);
        %!param(3) not relevant
        n    = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = 1/( sqrt(2*pi) * stdv ) * exp( -1/2 * ((x-mean)/stdv).^2 );
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 2  % Lognormal marginal distribution
        
        lambda = param(1); %mu
        zeta   = param(2); %sigma
        n      = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = 1./( sqrt(2*pi) * zeta * x ) .* exp( -1/2 * ((log(x)-lambda)/zeta).^2 );
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 3  % Gamma distribution
        
        lambda = param(1);
        k      = param(2);
        n      = param(4);
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = lambda * (lambda*x).^(k-1) / gamma(k) .* exp(-lambda*x);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 4  % Shifted exponential distribution
        
        lambda = param(1);
        x_zero = param(2);
        n      = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = lambda * exp(-lambda*(x-x_zero));
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 5  % Shifted Rayleigh distribution
        
        a      = param(1);
        x_zero = param(2);
        n      = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = (x-x_zero)/a^2 .* exp(-0.5*((x-x_zero)/a).^2);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 6  % Uniform distribution
        
        a = param(1);
        b = param(2);
        n = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = 1 / (b-a);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 7  % Beta distribution - NO powered option is avaialable!
        
        q = param(1);
        r = param(2);
        a = param(3);
        b = param(4);
        
        p = (x-a).^(q-1) .* (b-x).^(r-1) / ( (gamma(q)*gamma(r)/gamma(q+r)) * (b-a)^(q+r-1) );
        
    case 8  % Chi-square distribution
        
        nu     = param(1);
        lambda = 0.5;
        k      = nu/2;
        n      = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = lambda * (lambda*x).^(k-1) .* exp(-lambda*x) / gamma(k);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 11  % Type I largest value distribution ( same as Gumbel distribution )
        
        u_n = param(1);
        a_n = param(2);
        %!param(3) not relevant
        n   = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = a_n * exp( -a_n*(x-u_n) - exp(-a_n*(x-u_n)) );
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 12  % Type I smallest value distribution
        
        u_1 = param(1);
        a_1 = param(2);
        %!param(3) not relevant
        n   = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = a_1 * exp( a_1*(x-u_1) - exp(a_1*(x-u_1)) );
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 13  % Type II largest value distribution
        
        u_n = param(1);
        k   = param(2);
        %!param(3) not relevant
        n   = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = k/u_n * (u_n./x).^(k+1) .* exp(-(u_n./x).^k);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 14  % Type III smallest value distribution
        
        u_1     = param(1);
        k       = param(2);
        epsilon = param(3);
        n       = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = k/(u_1-epsilon) * ((x-epsilon)/(u_1-epsilon)).^(k-1) ...
                .* exp(-((x-epsilon)/(u_1-epsilon)).^k);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 15  % Gumbel distribution ( same as type I largest value distribution )
        
        u_n = param(1);
        a_n = param(2);
        %!param(3) not relevant
        n   = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = a_n * exp( -a_n*(x-u_n) - exp(-a_n*(x-u_n)) );
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 16  % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
        
        u_1 = param(1);
        k   = param(2);
        %!param(3) not relevant
        n   = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = k/u_1 * (x/u_1).^(k-1) .* exp(-(x/u_1).^k);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 18  % (Reserved for Laplace marginal distribution)
        
    case 19  % (Reserved for Pareto marginal distribution)
        
    case 20  % Generalized extreme value (GEV) distribution
        k       = param(1);
        sigma   = param(2);
        mu      = param(3);
        n       = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = gevpdf(x, k, sigma, mu);
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end
        
    case 25 % Three-parameter lognormal (LN3) distribution
        shape   = param(1);
        scale   = param(2);
        thres   = param(3);
        n       = param(4);
        
        if n == 1 % not powered distribution, to avoid finite difference calculation
            p = lognorm3pdf(x, shape, scale, thres, 'par');
        else % powered distribution
            Fn = @(x) ferum_cdf(type, x, param);
            
            p = cfd(Fn, x);
        end    
        
    case 30 % sample based custom distibution - typically for non-parametric distributions
        % powered option is valid, its effect is directly (numerically) incorporated into the kernel distribution (kernel_sample_*.mat)
        ID      = param(1);
        %!param(2) not relevant
        %!param(3) not relevant
        %!param(4) not relevant
        
        p       = custom_pdf(x, ID, 'sample');
        
    case 31 % vector based custom distibution - typically for non-parametric distributions
        % powered option is not valid
        ID      = param(1);
        %!param(2) not relevant
        %!param(3) not relevant
        %!param(4) not relevant
        
        p       = custom_pdf(x, ID, 'point');
        
    case 51  % Truncated normal marginal distribution; NO powered option is avaialable!
        
        mean = param(1);
        stdv = param(2);
        xmin = param(3);
        xmax = param(4);
        
        p = 1/(normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) * 1/( sqrt(2*pi) * stdv ) * exp( -1/2 * ((x-mean)/stdv).^2 );
        Imin = find(x<xmin); p(Imin) = 0;
        Imax = find(x>xmax); p(Imax) = 0;
        
    otherwise
        
end

end
