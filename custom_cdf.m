% Cumulative distrtibution function of a random variable
%
% P = CUSTOM_CDF(x, ID, type)
%
%INPUT
% x      - value at which we are interested in the cdf value, P(X<x)
% ID     - ID value of the distribution (specified in probdata.marg)
% type   - definition type of the custom distribution, 'sample' or 'point'
%
%OUTPUT
% P      - probability of x not being exceeded
%
%NOTE:
%

function P = custom_cdf(x, ID, type)

persistent S

switch (lower(type))
    
    case 'sample'
%         load(['tmp\kernel_sample_', num2str(ID), '.mat'], 'x_grid', 'cdf')
%         P = interp1(x_grid, cdf, x);
%         %P = interp1(x_grid, cdf, x, 'pchip');
   case 'point'
        % distribution function is given by discrete points
        if isempty(S)
            load('tmp\vector_distr_struct.mat', 'S')
        end
        
        x_grid = S(ID).x_grid;
        cdf    = S(ID).cdf;
        
        P = interp1(x_grid, cdf, x);
    otherwise
        error(['Unknown type: ', type])
end

end





