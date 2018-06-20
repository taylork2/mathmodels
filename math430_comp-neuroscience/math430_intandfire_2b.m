%2014-01-23

% Integrate and fire model

clear all
close all

C = 1; %uF
El = -60; %mV
Gl = 0.1; %mS 
Iapp = 2; %nA 
Ek = -85; %mV

%integrate and fire model 
Vth = -50; %mV threshold 
Vrst = -65; %mV %reset 
dg_sra = 0.1; 
tau_sra = 100; %msec 

%time 
T = 1000;
dt = 0.1;
t = 0:dt:T;

% calculate tau = R*C
R = 1/Gl;
tau = R * C;

%initialize numerical solution 
V = zeros(1,length(t));
V(1) = El;

%initialize g_sra 
g_sra = zeros(1, length(t));
g_sra(1) = 0;

spkcnt = 0;
tspk = zeros(1);
ISI = zeros(1);

%numerically calculate 
for j=1:length(t)-1
    kv1 = ( -V(j) + El - R*g_sra(j)*( V(j)-Ek ) + R*Iapp ) / tau;
    kg1 = -g_sra(j)/tau_sra; 
    av = V(j)+kv1*dt;
    ag = g_sra(j)+kg1*dt;
    kv2 = (-Gl*(av-El)-ag*(av-Ek)+Iapp)/C;
    kg2 = -ag/tau_sra;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;
    g_sra(j+1) = g_sra(j) + (kg1+kg2)*dt/2;
    
    %calculate when spike occurs new g_sra 
    if V(j+1)>=Vth
       V(j)=60;
       V(j+1)=Vrst;
       g_sra(j+1) = g_sra(j) + dg_sra;  
       spkcnt = spkcnt+1;
       tspk(spkcnt)=t(j);
   end
end

for k=1:length(tspk)-1
    ISI(k) = tspk(k+1)-tspk(k);
end

%plotting the numerical soln
figure(1)
hold on 
subplot(2,1,1)
plot(t,V,'b','linewidth',2);
legend('Numerical V')
axis([0 T -80 80]);
set(gca,'fontsize',10);
ylabel('V[mV]');
title(['Integrate and fire model with additional current (tau_{sra} = ', num2str(tau_sra), ')'])

%graph g_sra 
subplot(2,1,2)
plot(t, g_sra, 'g');
set(gca, 'fontsize', 10);
xlabel('t[msec]')
ylabel('g_{sra}');
axis([0 T 0 0.1]);

