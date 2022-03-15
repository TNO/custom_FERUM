% Density function function for hardcoded in spearate m-files
%
% P = HARDCODED_PDF(x, ID, shift, scale)
%
%INPUT
% x         value at which we are interested in the cdf value, P(X<x)
% ID        ID value of the distribution (specified in probdata.marg)
% shift  shift parameter /default = 0/
% scale     scale parameter /default = 1/
%
%OUTPUT
% P         density value at x
%
%NOTE:

function [P, x_grid, pdf] = hardcoded_pdf(x, ID, shift, scale)

if nargin < 4
    scale = 1;
end

if nargin < 3
    shift = 0;
end

switch ID
    case 201 % wind
        x_grid      = wind_x_grid;
        pdf         = wind_pdf;
    case 301 % snow
        x_grid      = snow_x_grid;
        pdf         = snow_pdf;
    case 1001 % dummy1
        x_grid      = dummy1_x_grid;
        pdf         = dummy1_pdf;    
    otherwise
        error(['Unknown ID: ', num2str(ID)])
end

xx          = (x - shift)/scale;

P           = 1/scale*interp1(x_grid, pdf, xx);

end