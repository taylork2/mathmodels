h = 6.626*10^-34; %Plancks constant
e = 1.6*10^-19; %e-  
m = 9.109*10^-31; %e- mass
a = 0.1425*10^-9; 
L = 0.135;
V = 2400:50:5000; %voltage input
lamda = h./sqrt(2*m*e*V);
d1 = sqrt(3)/2*a;
d2 = 3/2*a;
n=1;
D1 = 2*L*sin(2*asin(n*lamda./(2*d1)));
D2 = 2*L*sin(2*asin(n*lamda./(2*d2)));


Vexp = [2.4, 2.9, 3.4, 3.9]*10^3;
D1exp = [3.4 2.9 2.7 2.4]*10^-2;
D2exp = [5.9 5.1 4.7 4.4]*10^-2;
hold on
plot(Vexp, D1exp, 'x:', V, D1, '-', V, D2, '--')
plot(Vexp, D2exp, 'x:', V, D1, '-', V, D2, '--')


%{
VIS = 1./sqrt(Vexp); %Inverse square root of Vexp 
coef1 = polyfit(VIS, D1exp, 1); %linear fitting 
s1 = coef1(1); %slope 
b1 = coef1(2); %intercept
Yfit1 = s1*VIS + b1;


coef2 = polyfit(VIS, D2exp, 1); %linear fitting 
s2 = coef2(1); %slope 
b2 = coef2(2); %intercept
Yfit2 = s2*VIS + b2; 

plot(VIS, D1exp, 'x', VIS, Yfit1, '--')
hold on 
plot(VIS, D2exp, 'x', VIS, Yfit2, '--')
hold off
d1exp = (2*L)/s1*(h/sqrt(2*m*e))
d2exp = (2*L)/s2*(h/sqrt(2*m*e))

err1 = abs(d1-d1exp)/d1 * 100
err2 = abs(d2-d2exp)/d2 * 100
%}
