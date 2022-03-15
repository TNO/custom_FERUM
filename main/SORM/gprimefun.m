function allg = gprimefun(ETA,Uprime1,R1,lsf,probdata,analysisopt,gfundata,femodel,randomfield)

% Evaluate the aggrigated limit-state function - Include rho - ||u||

global nfun

% Find number of random variables
marg = probdata.marg;
nrv = size(marg,1);

Uprime = Uprime1;
Uprime(nrv,:) = ETA;
allu = R1'*Uprime;

% Transform into original space
allx = u_to_x( allu , probdata );

% Evaluate limit-state function
[ allg, dummy ] = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);