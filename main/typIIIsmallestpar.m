function f = typIIIsmallestpar(x,meaneps,stdv)

f = (gamma(1+2/x)-(gamma(1+1/x))^2)^0.5 - (stdv/meaneps)*gamma(1+1/x);

% System of equations to be solved for the type III smallest distribution parameters
% f(1) = epsilon + (x(1)-epsilon)*gamma(1+1/x(2)) - mean;
% f(2) = (x(1)-epsilon)*(gamma(1+2/x(2))-(gamma(1+1/x(2)))^2)^0.5 - stdv;