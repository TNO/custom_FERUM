function  sormpfresults = sorm_pf(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

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

marg = probdata.marg;
nrv = size(marg,1);

if beta < 1
   k = 1 / abs(beta);
elseif beta < 3
   k = 1;
else
   k = 3 / abs(beta);
end

% Points of ordinates +beta
Uprime1 = [ [ -k*beta*eye(nrv-1) ; beta*ones(1,nrv-1) ] [ k*beta*eye(nrv-1) ; beta*ones(1,nrv-1) ] ];
ETA1 = Uprime1(nrv,:);

switch lower(gfundata(lsf).evaluator)
   case 'basic'
      tolx = 1e-5;
   otherwise
      tolx = 1e-3;
end

[ ETA, G ] = my_fzero_gprimefun_vectorized(ETA1,tolx,Uprime1,R1,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
ai = 2*(ETA-beta)/(k*beta)^2;

sormpfresults.uprimen_minus    = ETA([1:(nrv-1)]);
sormpfresults.G_minus          = G([1:(nrv-1)]);
sormpfresults.uprimen_plus     = ETA([nrv:2*(nrv-1)]);
sormpfresults.G_plus           = G([nrv:2*(nrv-1)]);



% Breitung modified by Hochenbichler / Rackwitz
kk = normpdf(beta) / normcdf(-beta);
pf2_breitung_m = normcdf(-beta) * prod( 1/2 * ( (1./(1+kk*ai([1:(nrv-1)])).^0.5) + (1./(1+kk*ai([nrv:2*(nrv-1)])).^0.5) ) );
if ~isreal(pf2_breitung_m)
    warning('The Breitung approximation (SORM) yielded to complex probabilities, the norm of that will be used in further analysis. The corresponding result might be disregarded!')
    pf2_breitung_m = norm(pf2_breitung_m);
end
    betag_breitung_m = -inv_norm_cdf(pf2_breitung_m);

sormpfresults.betag_breitung_m = betag_breitung_m;
sormpfresults.pf2_breitung_m   = pf2_breitung_m;

sormpfresults.uprimei_minus    = -k*beta*ones(1,(nrv-1));
sormpfresults.uprimei_plus     = k*beta*ones(1,(nrv-1));
sormpfresults.ai_minus         = ai([1:(nrv-1)]);
sormpfresults.ai_plus          = ai([nrv:2*(nrv-1)]);

sormpfresults.nfun             = nfun;
