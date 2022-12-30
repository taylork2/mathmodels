q1=1*10^-9; q2=-1*10^-9
esp0=8.85*10^-12;
k=1./(4*pi*esp0);
x=-2.1:0.2:2.1;y=-2.1:0.2:2.1;
[X,Y] = meshgrid(x,y);
r1=(X.^2 + (Y-a).^2).^0.5;
r2=(X.^2 + (Y+a).^2).^0.5;
E1 = k*q1./r1.^2; E2 = k*q2./r2.^2;
[px1] = cos(atan2(Y-a,X)); [py1] = sin(atan2(Y-a,X));
xcompl = E1.*px1; ycompl = E1.*py1;
[px2] = cos(atan2(Y+a,X)); [py2] = sin(atan2(Y+a,X));
xcomp2=E2.*px2; ycomp2=E2.*py2;
xcomp = xcomp1 + xcomp2; ycomp = ycomp1 + ycomp2;
quiver(X, Y, xcomp, ycomp, 1)
