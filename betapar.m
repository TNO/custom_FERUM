function f = betapar(q,a,b,mean,stdv)

r = (b-mean)/(mean-a)*q;
%f = ((b-a)/(q+r))*(q*r/(q+r+1))^0.5 - stdv;
f = abs( ((b-a)/(q+r))*(q*r/(q+r+1))^0.5 - stdv );          

% System of equations to be solved for the beta distribution parameters.
% f(1) = (a*r+b*q)/(r+q) - mean;
% f(2) = (b-a)/(r+q)*sqrt(r*q/(r+q+1)) - stdv;
