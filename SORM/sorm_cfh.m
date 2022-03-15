function  sormcfhresults = sorm_cfh(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

%     Finite Element Reliability Using Matlab, FERUM, Version 4.0, 2009 
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     A copy of the GNU General Public License is found in the file
%     <gpl-3.0.txt> following this collection of FERUM program files.
%     This license can be also found at: <http://www.gnu.org/licenses/>.    
%
%     For more information on FERUM, visit: <http://www.ifma.fr/FERUM/>


global nfun
nfun = 0;

beta   = analysisopt.formresults.beta;
dsptx  = analysisopt.formresults.dsptx;
dsptu  = analysisopt.formresults.dsptu;
alpha  = analysisopt.formresults.alpha;
G      = analysisopt.formresults.G;
grad_G = analysisopt.formresults.grad_G;

% Compute the Gram-Schmidt orthonormal matrix R1 if needed
if isfield(analysisopt.formresults,'R1')
   R1 = analysisopt.formresults.R1;
else
   R1 = orthonormal_matrix(analysisopt.formresults.alpha); analysisopt.formresults.R1 = R1;
end


hess_G = hessian(lsf,dsptx,dsptu,G,probdata,analysisopt,gfundata,femodel,randomfield);


marg = probdata.marg;
nrv = size(marg,1);

A = R1*hess_G*R1' / norm(grad_G);

kappa = eig(A([1:(nrv-1)],[1:(nrv-1)]));
% Looks for abberant complex values
I = find(abs(imag(kappa))<eps);
if ~isempty(I)
   kappa(I) = real(kappa(I));
end

% Breitung
pf2_breitung = normcdf(-beta) * prod(1./(1+beta*kappa).^0.5);
betag_breitung = -inv_norm_cdf(pf2_breitung);

% Breitung modified by Hochenbichler / Rackwitz
kk = normpdf(beta) / normcdf(-beta);

pf2_breitung_m = normcdf(-beta) * prod(1./(1+kk*kappa).^0.5);
betag_breitung_m = -inv_norm_cdf(pf2_breitung_m);

sormcfhresults.betag_breitung   = betag_breitung;
sormcfhresults.pf2_breitung     = pf2_breitung;

sormcfhresults.betag_breitung_m = betag_breitung_m;
sormcfhresults.pf2_breitung_m   = pf2_breitung_m;

sormcfhresults.kappa            = kappa;

sormcfhresults.nfun             = nfun;