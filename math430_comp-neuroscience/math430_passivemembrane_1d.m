%2017-09-18

% Passive membrane equation

clearvars;
close all

C = 1; %uF
El = -60; %mV
Gl = 0.1; %mS
I0 = 0.5; %nA
w = [ 1 5 10 20 100];
D = 0;

%time 
T = 1000; %msec
dt = 1; 
t = 0:dt:T;

% calculate tau = R*C
tau = 1/Gl * C;

%Numerical solution initialization
V = zeros(1,length(t));
V(1) = El;

%analytical solution initilization
A = I0;
freq = 2*pi*w(2)/1000;
k_an = A *(1/Gl) / (1 + freq^2*tau^2);
k_an2 = k_an * (sin(freq*t) - freq*tau*cos(freq*t));
V_an = El + k_an*freq*tau*exp(-t/tau) + k_an2;

figure(1)
hold on
%plot(t,V_an,'r','linewidth',1);

%calculate numerical 
for i=1:length(w)
    %assign new Iapp depending on w 
    Iapp = A * sin(2*pi*w(i)*t/1000);
    
    %runge kutta
    for j=1:length(t)-1
        kv1 = (-Gl*(V(j)-El)+Iapp(j))/C;
        av = V(j)+kv1*dt;
        kv2 = (-Gl*(av-El)+Iapp(j))/C;
        V(j+1) = V(j) + (kv1+kv2)*dt/2;  
    end
    
    plot(t,V,'linewidth',2);

end

%plotting extras 
legend('w=1', 'w=5', 'w=10', 'w=20', 'w=100');
set(gca,'fontsize',10);
xlabel('t[msec]')
ylabel('V[mV]');
title(['Passive membrane equation (Iapp=Asin(2*pi*w*t/1000))'])

% graph the input vs output freq
figure(2) 
plot(w, 2*pi*w/1000);
xlabel('input freq, w');
ylabel('output freq');
title('Input frequency vs. Output frequency')

%graph the input vs amp 
figure(3)
Amp = (1/Gl)*A ./ sqrt(1 + w.^2*tau^2);
plot(w, Amp);
xlabel('input frequency, w')
ylabel('output amplitude, A')
title('Input frequency vs. Amplitude')



