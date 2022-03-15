function allu = ssvm_sdu(nvar,nsample)

% s = sdu(nvar,nsample)
% SDU uniform distribution in the unit hypersphere
%
% Input:
%   nvar    : space dimension
%   nsample : nb. of samples
%
% Output:
%   s       : random sample (nvar,nsample)

allu = inv_norm_cdf(twister(nvar,nsample));
alla = allu ./ ( ones(nvar,1) * dot(allu,allu).^0.5 );

allu = twister(1,nsample);
allr = allu.^(1/nvar);

allu = ( ones(nvar,1) * allr ) .* alla;