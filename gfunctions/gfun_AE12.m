function g = gfun_AE12(x1,x2)

y1 = 3+0.1*(x1-x2).^2-(x1+x2)/sqrt(2);
y2 = 3+0.1*(x1-x2).^2+(x1+x2)/sqrt(2);
y3 = (x1-x2)+3.5*sqrt(2);
y4 = (x2-x1)+3.5*sqrt(2);

g = min([y1; y2; y3; y4],[],1);