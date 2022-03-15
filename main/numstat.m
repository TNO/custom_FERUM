% Calculate the mean and standard deviation of a distribution defined by probdata.marg
% numerical integration - adaptive quadrature - for non-special function, e.g. powered
%
%SYNOPSYS
% [mean, stdv] = NUMSTAT(margi, zmax)
%
function [mean, stdv] = numstat(margi, zmax)

pon = margi(11);
% integration limits
if margi(11) == 1;
    zmin = -zmax;
else
    zmin = norminv(ferum_cdf(margi(1), 0, margi(5:end))); %! WARNING
    %zmin = max(norminv(ferum_cdf(margi(1), 0, margi(5:end))), -zmax);
end
xmax = ferum_invcdf(margi(1), normcdf(zmax), margi(5:end));
xmin = ferum_invcdf(margi(1), normcdf(zmin), margi(5:end));

%area = integral(@(x) ferum_pdf(margi(1), x, margi(5:end)), xmin, xmax) + (1-pon);
mean = integral(@(x) ferum_pdf(margi(1), x, margi(5:end)).*x, xmin, xmax);
stdv = sqrt(integral(@(x) ferum_pdf(margi(1), x, margi(5:end)).*(x-mean).^2, xmin, xmax) + (1-pon)*mean^2);

end