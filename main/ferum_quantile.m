function xq = ferum_quantile(x,p)

% x  : row vector
% p  : scalar value ( 0 <= p <= 1 )
% xq : p-th quantile

x = sort(x,2);
n = length(x);

q = [0 100*(0.5:(n-0.5))/n 100];
xx = [x(1) x(1:n) x(n)];
xq = interp1q(q',xx',100*p);
