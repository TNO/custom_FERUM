function G = gprime(lsf,Uprime,R1,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun

marg = probdata.marg;
nrv = size(marg,1);

U = R1'*Uprime;

X = zeros(nrv,size(Uprime,2));
G = zeros(1,size(Uprime,2));

X = u_to_x(U,probdata);

[ G, dummy ] = gfun(lsf,X,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
