% Inverse cumulative distrtibution function of a custom random variable
%
% x = CUSTOM_INVCDF(u, ID, type)
%
%INPUT
% P      - probability of x not being exceeded, P(X<x)
% ID     - ID value of the distribution (specified in probdata.marg)
% type   - definition type of the custom distribution, 'sample' or 'point'
%
%OUTPUT
% x      - value at which the probability on non-exceedance equals to P
%
%NOTE:
%

function x = custom_invcdf(u, ID, type)

persistent S

switch (lower(type))
    
    case 'sample'
        load(['tmp\kernel_sample_', num2str(ID), '.mat'], 'x_grid', 'cdf')
        x = interp1(cdf, x_grid, normcdf(u));
        %x = interp1(cdf, x_grid, P, 'pchip');
    case 'point'
        % distribution function is given by discrete points
        % distribution function is given by discrete points
        if isempty(S)
            load('tmp\vector_distr_struct.mat', 'S')
        end
        
        x_grid = S(ID).x_grid;
        cdf    = S(ID).cdf;
        
% %         nu = normcdf(u);
% %         if nu == 1 || nu == 0
% %             if exist('mp.m', 'file') == 2
% %                 mp.Digits(68);
% %                 u = mp(u);
% %                 x = interp1(cdf, x_grid, normcdf(u));
% %                 x = double(x);
% %             else
% %                 x = interp1(cdf, x_grid, normcdf(u));
% %             end
% %         else    
% %             x = interp1(cdf, x_grid, normcdf(u));
% %         end
        
        x = interp1(cdf, x_grid, normcdf(u));

    otherwise
        error(['Unknown type:', type])
end

end





