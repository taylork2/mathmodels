lambda = [447,486,502,546,656]*10^-9; %wavelengths
c = 3 * 10^8; %speed of light
e = 1.6 * 10^-19; %electron charge
h = 6.63*10^-34; %planck constant
phi = 2.7; 
frequenciesExp = c./lambda; %experimental frequencies
Xexp = frequenciesExp *10^19*10^-34;
WF = phi * 10^-19/e; %work function
vRead = [1.075 1.048 1.044 1.129 1.069;
    0.505 0.441 0.476 0.454 0.455;
    0.602 0.616 0.623 0.588 0.616;
    0.630 0.582 0.580 0.587 0.601;
    0.314 0.359 0.314 0.365 0.301]; %read from multimeter
vExp = mean(vRead'); %take mean for accuracy
K = e*vExp;
Yexp = K*10^19;
coeff = polyfit(Xexp,Yexp,1);
m = coeff(1) %slope
b = coeff(2) %intercept
Ybestfit = m*Xexp + b;
Xth = 0.1:0.001:0.8;
Yth = (h/10^-34).*Xexp - phi;
Err = abs(m*10^-34-h)/h*100
Yth_exp = (h/10^-34).*Xexp - phi;
Err_rel = (Yexp - Yth_exp)./Yth_exp*100


