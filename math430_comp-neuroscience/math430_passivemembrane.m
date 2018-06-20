%1b

% Passive membrane equation
clearvars;
close all

C = 1; %uF
El = -60; %mV
Gl = 0.1; %mS
Iapp = -0.5; %nA
D = 0; 

T = 400; %end time [msec]
dt = 0.1; 
t = 0:dt:T;

% calculate tau = R*C
tau = 1/Gl * C;

%initialize numerical solution 
V = zeros(1,length(t));
V(1) = El;

%runge kutta - calc numerical soln 
for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp)/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp)/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;   
end

%analytical solution 
Vinf = 1/Gl*Iapp;
V_an = Vinf+El-Vinf*exp(-t/tau);

%plot numerical & analytical solution 
figure(1)
hold on
plot(t,V,'b','linewidth',2);
plot(t,V_an,'r','linewidth',1);
legend('Numerical V', 'Analytical V');
axis([0 T -80 -30]);
set(gca,'fontsize',10);
xlabel('t[msec]')
ylabel('V[mV]');
title(['Passive membrane equation (Iapp=', num2str(Iapp), ')'])

%absolute value error
V_dif = abs(V-V_an);

%plot error 
figure(2)
plot(t,V_dif,'g','linewidth',1);
title(['Error Iapp=',num2str(Iapp)]);
axis([0 T 0 1]);
xlabel('t[msec]');
ylabel('Error[mV]');


