% Two-point central fintite difference method to approximate derivative
%
%SYNOPSYS
% p = CFD(Fn, x)
%
% Used to approximate the pdf at a given point from the cdf using finite difference method (central difference with two points)
% [Karen A. Kopecky (2007). Lecture Notes. http://www.karenkopecky.net/Teaching/eco613614/Notes_NumericalDifferentiation.pdf]
% ~10e-11 accuracy (O(h^2))
%
%INPUT
% Fn - handler of the cdf function
% x  - point(s) where we are interested in the value of pdf corresponding to Fn cdf, /can be vector/
%
%OUTPUT
% p  - numerically estimated pdf function at point x /same size and type as x/
%
%See also
% oimatrix

function p = cfd(Fn, x, h)
    if nargin < 3
        h   = (eps)^(1/3)*arrayfun(@(y) max([y, 1]), abs(x));
    end
    xph = x + h;
    xmh = x - h;
    dx  = xph - xmh; % to account for rounding error

    p = (Fn(xph) - Fn(xmh))./(dx);
end

