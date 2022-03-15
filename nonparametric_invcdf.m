% Inverse cumulative distrtibution function for nonparametric distribution (defined
% with discrete points!)
%
% x = NONPARAMETRIC_INVCDF(P, ID, shift, scale)
%
%INPUT
% x         value at which we are interested in the cdf value, P(X<x)
% ID        ID value of the distribution (specified in probdata.marg)
% shift  shift parameter /default = 0/
% scale     scale parameter /default = 1/
%
%OUTPUT
% x         x at P
%
%NOTE:

function x = nonparametric_invcdf(u, ID, shift, scale)


if nargin < 4
    scale = 1;
end

if nargin < 3
    shift = 0;
end

persistent S

% distribution function is given by discrete points
if isempty(S)
    load('tmp\vector_distr_struct.mat', 'S')
end

x_grid      = S(ID).x_grid;
cdf         = S(ID).cdf;

x           = scale*interp1(cdf, x_grid, normcdf(u)) + shift;

end