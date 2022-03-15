% Density function for nonparametric distribution (defined
% with discrete points!)
%
% P = NONPARAMETRIC_PDF(x, ID, location, scale)
%
%INPUT
% u         (due to numerical issues)
% ID        ID value of the distribution (specified in probdata.marg)
% location  location parameter /default = 0/
% scale     scale parameter /default = 1/
%
%OUTPUT
% P         density value at x
%
%NOTE:

function P = nonparametric_pdf(x, ID, shift, scale)

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
pdf         = S(ID).pdf;

xx          = (x - shift)/scale;

P           = 1/scale*interp1(x_grid, pdf, xx);

end