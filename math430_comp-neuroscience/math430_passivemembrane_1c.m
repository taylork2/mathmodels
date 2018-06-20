%2017-09-18

% Passive membrane equation

clearvars;
close all

C = 1; %uF
El = -60; %mV
Gl = 0.1; %mS
I0 = 0.5; %nA
D = 0;

%time 
T = 400; %msec
dt = 0.1; 
t = 0:dt:T;

% calculate tau = R*C
tau = 1/Gl * C;

%Heaviside equation, used to initialize Iapp
ti = 100;
tf = 200;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

%Numerical solution initialization
V = zeros(1,length(t));
V(1) = El;

%analytical solution initilization
Vinf = 1/Gl*I0;
V_an = V;

%Initializing Iapp values 
Iapp = zeros(1, length(t));
for i=1:length(t)
    Iapp(i) = I0*H(i); %nA
end

%runge kutta
for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp(j)*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp(j)*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;  
    
    %analytical
    V_an(j+1) = El + Vinf * ( heaviside(t(j)-ti) * (1 - exp(-(t(j)-ti)/tau)) - heaviside(t(j)-tf) * (1 - exp(-(t(j)-tf)/tau)) );

end

%absolute value 
V_dif = abs(V-V_an);

%plotting numerical & analytical soln 
figure(1)
hold on
plot(t,V,'b','linewidth',2);
plot(t,V_an,'r','linewidth',1);
legend('Numerical V', 'Analytical V');
axis([0 T -80 -20]);
set(gca,'fontsize',10);
xlabel('t[msec]')
ylabel('V[mV]');
title(['Passive membrane equation (Iapp=', num2str(I0), ')'])

%plotting error 
figure(2)
plot(t,V_dif,'g','linewidth',1);
title(['Error Iapp=', num2str(I0)])
xlabel('t[msec]')
ylabel('Error[mV]');
