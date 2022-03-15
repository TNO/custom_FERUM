function [mn, mx, d] = plot_histo(i,xi,N,odd,scale)

mn = min(xi);
mx = max(xi);
d = (mx - mn)/N*2;
e = floor(log(d)/log(10));
m = floor(d/10^e);
if m > 5
   m = 5;
elseif m > 2
   m = 2;
end
d = m * 10^e;
mn = (floor(mn/d)-1)*d - odd*d/2;
mx = (ceil(mx/d)+1)*d + odd*d/2;
limits = mn:d:mx;

f = zeros(1,length(limits)-1);
for j = 1:length(limits)-1
   f(j) = sum( xi>=limits(j) & xi<limits(j+1) );
end

xx = [limits; limits; limits];
xx = xx(:);
xx = xx(2:length(xx)-1);
yy = [f*0; f; f];
yy = [yy(:); 0];
if scale, yy = yy/length(xi)/d; end

figure(i)
clf
plot(xx,yy,'b-')
hold on
plot(limits,limits*0)
set(gca,'FontSize',18);