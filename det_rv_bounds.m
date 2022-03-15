% Determine lower and upper bounds of random variables
%
%SYNOPSYS
% probdata = DET_RV_BOUNDS(probdata)
%
%

function probdata = det_rv_bounds(probdata, umax)

if nargin < 2
   umax = []; 
end

marg    = probdata.marg;
n_rv    = size(marg,1);
bounds  = nan(n_rv, 2);

for ii = 1:n_rv
    dist_ii = marg(ii,1);
    switch dist_ii
        case 1  % Normal distribution
            bounds_ii = [-Inf, Inf];
        case 2  % Lognormal distribution
            bounds_ii = [0, Inf];
        case 3  % Gamma distribution
            
        case 4  % Shifted exponential distribution
            
        case 5  % Shifted Rayleigh distribution
            
        case 6  % Uniform distribution
            a           = marg(ii,5);
            b           = marg(ii,6);
            bounds_ii    = [a, b];
        case 7  % Beta distribution
            
        case 8 % Chi-square distribution
            
        case 11 % Type I largest value distribution ( same as Gumbel distribution )
            bounds_ii = [-Inf, Inf];
        case 12 % Type I smallest value distribution
            
        case 13 % Type II largest value distribution
            
        case 14 % Type III smallest value distribution
            
        case 15 % Gumbel distribution ( same as type I largest value distribution )
            bounds_ii = [-Inf, Inf];
        case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
            
        case 18 % (Reserved for Laplace distribution)
        case 19 % (Reserved for Pareto distribution)
            
        case 20 % Generalized extreme value (GEV) distribution
            
        case 30 % sample based custom distribution
            
        case 31 % vector based custom distribution
            ID          = marg(ii,5);
            load(['tmp\vector_distr_',num2str(ID),'.mat'], 'x_grid')
            bounds_ii    = [min(x_grid), max(x_grid)];
        case 51 % Truncated normal marginal distribution
            
        otherwise
            error(['Unknown distribution type: ', dist_ii])
    end
    bounds(ii,:)    = bounds_ii;
    bounds_ii       = [];
end

probdata.bounds = bounds;

end