function ally = sigmafun(allr,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield)

% Evaluate the aggrigated limit-state function - Include rho - ||u||

global nfun

% Find number of random variables
nrv = size(alla,1);

% Transform into original space
allx = u_to_x( ( ones(nrv,1)*allr ) .* alla , probdata );

% Evaluate limit-state function
[ allg, dummy ] = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);

ally = min( allg , rho - allr );