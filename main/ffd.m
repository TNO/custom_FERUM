% Two-point forward fintite difference method to approximate derivative
%
%SYNOPSYS
% grad = FFD(F, x)
%
% [Karen A. Kopecky (2007). Lecture Notes. http://www.karenkopecky.net/Teaching/eco613614/Notes_NumericalDifferentiation.pdf]
% ~10e-8 accuracy (O(h))
%
%INPUT
% F  - handler of the function /vector valued/
% x  - point where we are interested in the gradient
%
%OUTPUT
% grad  - gradient vector - finite forward approximate derivative at point x /same size and type as x/
%
%NOTE
% Works only for one point at once (not vectorized)!
%
%See also
% cfd

function grad = ffd(F, x, h)

x = x(:).';

if nargin < 3
    h   = (eps)^(1/2)*arrayfun(@(y) max([y, 1]), abs(x));
end

nvar = length(x);
% evaluate at x
F0 = F(x);

xx  = repmat(x,nvar,1);
xph = xx + diag(h);
dx  = diag(xph - xx); % to account for rounding error

grad = nan(size(x));
% who knows if F is vectorized or not
for ii = 1:nvar
    grad(ii)   = (F(xph(ii,:)) - F0)/(dx(ii));
end
% disp('ffd')
end

