syms o theta M hstar n m h x1 x2 x3 x4 x5 s v v2 gamma dphi1 dphi2

dphi = @(h) gamma * (-n*hstar^n/h^(n+1) + m*hstar^m/h^(m+1));

g1 = ((x4+x3)^3) *(x4-x3) - ((x3+x2)^3)*(x3-x2);

g2 = ((x4+x3)^3) * (x5-3*x4+3*x3-x2) - ((x3+x2)^3) * (x4-3*x3+3*x2-x1);

g3 = (x3+x4)^3*(x4-x3)*dphi1 - (x2+x3)^3*(x3-x2)*dphi2;

F = x3-(s*g1 + v*g2 + v2*g3);

diff(F, x1);
diff(F, x2)
diff(F, x3)
diff(F, x4)
diff(F, x5);

% diff(F,x3)
% 
% diff(F,x2)
% 
% diff(F,x4)