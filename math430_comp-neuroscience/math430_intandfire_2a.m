%2014-01-23

% Integrate and fire model

clear all
close all

C = 1; %uF
El = -60; %mV
Gl = 0.1; %mS 
Iapp = 2; %nA %plug in 1, 1.01, and 2

%integrate and fire model 
Vth = -50; %mV threshold 
Vrst = -65; %mV %reset 

%time 
T = 1000;
dt = 0.1;
t = 0:dt:T;

% calculate tau = R*C
tau = 1/Gl * C;

%Heaviside equation 
ti = dt;
tf = 10000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

%initialize numerical solution 
V = zeros(1,length(t));
V(1) = El;

%calculate r_isi
Vinf = 1/Gl*Iapp;
t_isi = tau * log( (Vrst - Vinf - El) / (Vth - Vinf - El) );
r_isi = 1/t_isi;
r_isi

spkcnt = 0;
tspk = zeros(1);
ISI = zeros(1);

%numerically calculate 
for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;
    if V(j+1)>=Vth
       V(j)=60;
       V(j+1)=Vrst;
       spkcnt = spkcnt+1;
       tspk(spkcnt)=t(j);
   end
end

for k=1:length(tspk)-1
    ISI(k) = tspk(k+1)-tspk(k);
end

%analytical solution 
V_an = Vinf+El+(V(1)-El-Vinf)*exp(-t/tau);

%plotting the analytical and numerical soln
figure(1)
hold on
plot(t,V,'b','linewidth',2);
plot(t,V_an,'r','linewidth',1);
legend('Numerical V', 'Analytical V')
axis([0 T -80 80]);
set(gca,'fontsize',10);
xlabel('t[msec]')
ylabel('V[mV]');
title(['Integrate and fire model (Iapp = ', num2str(Iapp), ')'])

%absolute value 
V_dif = abs(V-V_an);

%plot error 
figure(2)
plot(t,V_dif,'g','linewidth',1);
title(['Error Iapp=',num2str(Iapp)]);
axis([0 T 0 80]);
xlabel('t[msec]');
ylabel('Error[mV]');


