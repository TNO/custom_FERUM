function f = typIIlargestpar(x,mean,stdv)

f = (gamma(1-2/x)-(gamma(1-1/x))^2)^0.5 - (stdv/mean)*gamma(1-1/x);

% System of equations to be solved for the type II largest distribution parameters
% f(1) = x(1)*gamma(1-1/x(2)) - mean;
% f(2) = x(1)*(gamma(1-2/x(2))-(gamma(1-1/x(2)))^2)^0.5 - stdv;

