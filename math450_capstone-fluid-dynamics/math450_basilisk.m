clear all;
close all; 

fid=fopen('math450_wave128.txt');
s=textscan(fid,'%f %f');
fclose(fid);
t=s{1};
h=s{2};

E=0.01; 
h0=.1;

y=log(abs((0.9+h)/(E*h0))); %add 0.9 to shift

coef = polyfit(t, y, 1);
slope = coef(1)

figure(1)
plot(t,y);
title('Basilisk Linear Stability at x_{max} vs t');
xlabel('t');
ylabel('ln|(h(x_{max},t)-h0)/E*h0)|');

figure(2)
plot(t,h);
title('Basilisk x_{max} vs t');
xlabel('t');
ylabel('x_{max}');