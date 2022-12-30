%plot mass vs T^2 and regression line 
Texp = [11.3, 11.96, 12.60, 13.19, 13.72, 14.28]/20; %Period
Mexp = [0.025,0.035,0.045, 0.055, 0.065, 0.075]; %Additional mass
mm=0.01:0.001:0.1;
Tsquare = Texp.^2;
plot(Mexp, Tsquare, 'o')
coeff = polyfit(Mexp, Tsquare,1)
a = coeff (1);
b = coeff (2);
MM = 0.01:0.001:0.1;
ybestfit = a*MM +b;
plot(Mexp, Tsquare, 'o', mm, ybestfit, ':')
xlabel('Mass (kg)')
ylabel('Tsquare')
title('Tsquare vs Mass')
legend('Tsquare', 'ybestfit')
