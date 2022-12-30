%lab lab s
z = 1;
k = 9*10^9;
e = 1.6*10^-19; %fundamental charge
me = 9.1*10^-31; %mass of electron
h = 6.626*10^-34; %planck constant
n = 1:7; %orbits
m1 = 1;
M1 = 1;
mu_hydrogen = m1*M1/(m1+M1);
m2 = 1;
M2 = 2;
mu_deuterium = m2*M2/(m2+M2);
Diff1 = (mu_hydrogen-1)/1*1;
Diff2 = (mu_deuterium-1)/1*1;
%wavelength of stationary electron mass

E = (4*pi^2*me*k^2*e^4*z^2)./(2*h^2*n.^2);
plot(n,E,'*'),xlabel 'n', ylabel 'E'
axis([0 8 0 1*10^-18])
ni = [2,3,4,5,6,3,4,5,6,4,5,6];
nf = [1,1,1,1,1,2,2,2,2,3,3,3];
Ei = (4*pi^2*me*k^2*e^4*z^2)./(2*h^2*ni.^2);
Ef = (4*pi^2*me*k^2*e^4*z^2)./(2*h^2*nf.^2);
deltaE = (Ef-Ei)./e;
c = 3*10^8; %speed of light
lamda = (h*c)./(Ef-Ei);
lamda' %display lamda'

%wavelength based on the reduced mass

mu = me/(1+me/M1);
Ei_mu = (4*pi^2*me*k^2*e^4*z^2)./(2*h^2*ni.^2);
Ef_mu = (4*pi^2*me*k^2*e^4*z^2)./(2*h^2*nf.^2);
En = (4*pi^2*me*k^2*e^4*z^2)./(2*h^2*n.^2);
plot(n,En), xlabel 'n', ylabel 'En'

deltaE_mu = Ei_mu - Ef_mu;
lambda_mu = (h*c)./deltaE_mu;
lambda_mu'