function u = x_to_u(x,probdata)
%X_TO_U    Transformation between 'original' and 'standard normal' space
%
%   u = x_to_u(x,probdata)
%
%   Perform the transformation between original space and standard
%   normal space for a random variable with given probability distribution.
%
%   Output: - u           = values of r.v.'s in standard normal space
%   Input:  - x           = values of r.v.'s in original space
%           - probdata

%! NOTE!
%! the same procedure is done in every cases, this whole case stuff could be replaced

transf_type = probdata.transf_type;
marg = probdata.marg;
nrv = size(marg,1);

switch transf_type
    
    case 1
        
    case 2
        
        %u = zeros(nrv,1);
        %for i = 1 : nrv
        %    switch marg(i,1)
        %       case 1  % Normal marginal distribution
        %          u(i) =  ( x(i) - marg(i,2) ) / marg(i,3) ;
        %       case 2  % Lognormal marginal distribution
        %          zeta = sqrt( log( ( 1 + (marg(i,3)/marg(i,2))^2 ) ) );
        %          lambda = log( marg(i,2) ) - 0.5 * zeta^2;
        %          u(i) =  ( log(x(i)) - lambda ) / zeta ;
        %       case 4  % Shifted exponential marginal distribution
        %          u(i) = inv_norm_cdf( ferum_cdf(4,x(i),marg(i,2),marg(i,3)) );
        %       case 5  % Shifted Rayleigh marginal distribution
        %          u(i) = inv_norm_cdf( ferum_cdf(5,x(i),marg(i,2),marg(i,3)) );
        %       case 6  % Uniform marginal distribution
        %          u(i) = inv_norm_cdf( ferum_cdf(6,x(i),marg(i,2),marg(i,3)) );
        %       case 11 % Type I Largest Value or Gumbel marginal distribution
        %          u(i) = inv_norm_cdf( ferum_cdf(11,x(i),marg(i,2),marg(i,3)) );
        %       otherwise
        %    end
        %end
        
    case 3
        
        iLo = probdata.iLo;
        for i = 1 : nrv
            switch marg(i,1)
                case 1  % Normal marginal distribution
                    mean = marg(i,2);
                    stdv = marg(i,3);
                    n       = marg(i,8);
                    if n == 1
                        z(i,:) = (x(i,:) - mean)/stdv ;
                    else
                        z(i,:) = inv_norm_cdf( ferum_cdf(1, x(i,:), marg(i,5:8)) );
                    end
                case 2  % Lognormal marginal distribution
                    lambda  = marg(i,5);
                    zeta    = marg(i,6);
                    n       = marg(i,8);
                    mu      = lambda;
                    sigma   = zeta;
                    if n == 1
                        z(i,:) = (log(x(i,:)) - mu)/sigma;
                    else
                        z(i,:) = inv_norm_cdf( ferum_cdf(2, x(i,:), marg(i,5:8)) );
                    end
                    %z(i,:) =  ( log(x(i,:)) - lambda ) / zeta ;
                case 3  % Gamma distribution
                    %lambda = marg(i,5);
                    %k = marg(i,6);
                    %z(i) = inv_norm_cdf( gammainc(lambda*x(i),k) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(3,x(i,:),marg(i,5:8)) );
                case 4  % Shifted exponential distribution
                    %lambda = marg(i,5);
                    %x_zero = marg(i,6);
                    %z(i) = inv_norm_cdf( 1 - exp( -lambda*(x(i)-x_zero) ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(4,x(i,:),marg(i,5:8)) );
                case 5  % Shifted Rayleigh distribution
                    %a = marg(i,5);
                    %x_zero = marg(i,6);
                    %z(i) = inv_norm_cdf( 1 - exp( -0.5*( (x(i)-x_zero)/a )^2 ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(5,x(i,:),marg(i,5:8)) );
                case 6  % Uniform distribution
                    %a = marg(i,5);
                    %b = marg(i,6);
                    %z(i) = inv_norm_cdf( (x(i)-a) / (b-a) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(6,x(i,:),marg(i,5:8)) );
                case 7  % Beta distribution
                    %q = marg(i,5);
                    %r = marg(i,6);
                    %a = marg(i,7);
                    %b = marg(i,8);
                    %% reduce x to the interval [0,1]
                    %x01 = (x(i)-a) / (b-a);
                    %z(i) = inv_norm_cdf( betainc(x01,q,r) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(7,x(i,:),marg(i,5:8)) );
                case 8 % Chi-square distribution
                    %lambda = 0.5;
                    %nu = marg(i,5);
                    %k = nu/2 ;
                    %z(i) = inv_norm_cdf( gammainc(lambda*x(i),k) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(8,x(i,:),marg(i,5:8)) );
                case 11 % Type I largest value distribution ( same as Gumbel distribution )
                    %u_n = marg(i,5);
                    %a_n = marg(i,6);
                    %z(i) = inv_norm_cdf( exp( -exp( -a_n*(x(i)-u_n) ) ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(11, x(i,:), marg(i,5:8)) );
                case 12 % Type I smallest value distribution
                    %u_1 = marg(i,5);
                    %a_1 = marg(i,6);
                    %z(i) = inv_norm_cdf( 1 - exp( -exp( a_1*(x(i)-u_1) ) ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(12,x(i,:),marg(i,5:8)) );
                case 13 % Type II largest value distribution
                    %u_n = marg(i,5);
                    %k = marg(i,6);
                    %z(i) = inv_norm_cdf( exp( -(u_n/x(i))^k ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(13,x(i,:),marg(i,5:8)) );
                case 14 % Type III smallest value distribution
                    %u_1 = marg(i,5);
                    %k = marg(i,6);
                    %epsilon = marg(i,7);
                    %z(i) = inv_norm_cdf( 1 - exp( -( (x(i)-epsilon) / (u_1-epsilon) )^k ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(14,x(i,:),marg(i,5:8)) );
                case 15 % Gumbel distribution ( same as type I largest value distribution )
                    %u_n = marg(i,5);
                    %a_n = marg(i,6);
                    %z(i) = inv_norm_cdf( exp( -exp( -a_n*(x(i)-u_n) ) ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(15,x(i,:),marg(i,5:8)) );
                case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
                    %u_1 = marg(i,5);
                    %k = marg(i,6);
                    %z(i) = inv_norm_cdf( 1 - exp( -(x(i)/u_1)^k ) );
                    z(i,:) = inv_norm_cdf( ferum_cdf(16,x(i,:),marg(i,5:8)) );
                case 18 % (Reserved for Laplace distribution)
                case 19 % (Reserved for Pareto distribution)
                    
                case 20 % Generalized extreme value (GEV) distribution
                    z(i,:) = inv_norm_cdf( ferum_cdf(20,x(i,:),marg(i,5:8)) );
                case 25 % Three-parameter lognormal (LN3) distribution
                    z(i,:) = inv_norm_cdf( ferum_cdf(25,x(i,:),marg(i,5:8)) );
                case 30 % sample based custom distibution
                    z(i,:) = inv_norm_cdf( ferum_cdf(30,x(i,:),marg(i,5:8) ));
                case 31 % sample based custom distibution
                    z(i,:) = inv_norm_cdf( ferum_cdf(31,x(i,:),marg(i,5:8) ));
                case 32 % sample based custom distibution
                    z(i,:) = inv_norm_cdf( ferum_cdf(32,x(i,:),marg(i,5:8) ));
                case 33 % sample based custom distibution
                    z(i,:) = inv_norm_cdf( ferum_cdf(33,x(i,:),marg(i,5:8) ));    
                case 51 % Truncated normal marginal distribution
                    z(i,:) = inv_norm_cdf( ferum_cdf(51,x(i,:),marg(i,5:8)) );
                otherwise
            end
        end
        u = iLo * z;
        
    case 4
        
        %u = zeros(nrv,1);
        %for i = 1:nrv
        %   u(i) = inv_norm_cdf( cdf_morgenstern(x,marg,i) );
        %end
        
    otherwise
        
end

