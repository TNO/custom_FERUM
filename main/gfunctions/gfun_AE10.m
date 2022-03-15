function g = gfun_AE10(x1,x2,x3)

g = zeros(1,length(x3));
[dummy,J] = find( x3 <= 5 );
g(J) = x1(J)-x2(J)-x3(J);
[dummy,J] = find( x3 > 5 );
g(J) = x3(J)-x2(J);