function [ dbeta_dthetag, dpf_dthetag ] = sensitivities_wrt_thetag(lsf,beta,x,G,grad_G,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun

if ~isfield(analysisopt,'dthetagpara')
   switch lower(gfundata(lsf).evaluator)
      case 'basic'
         analysisopt.dthetagpara = 100000;
      otherwise
         analysisopt.dthetagpara = 100;
   end
end

thetag = gfundata(lsf).thetag;
nthetagv = length(thetag);

all_thetag = zeros(nthetagv,nthetagv);
allh = zeros(1,nthetagv);

original_thetag = thetag;
for j = 1:nthetagv
   thetag = original_thetag;
   if thetag(j) == 0
      allh(j) = 1/analysisopt.dthetagpara;
   else
      allh(j) = thetag(j)/analysisopt.dthetagpara;
   end
   thetag(j) = thetag(j) + allh(j);
   all_thetag(:,j) = thetag;
end

gfundata(lsf).thetag = all_thetag;

[ allG, dummy ] = gfun(lsf,x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);

% beta sensitivities w.r.t. limit-state parameters
dbeta_dthetag = ( ( allG(1:nthetagv) - G ) ./ allh(1:nthetagv) )' / norm(grad_G);

% pf sensitivities w.r.t. limit-state parameters
dpf_dthetag = -normpdf(beta) * dbeta_dthetag;
