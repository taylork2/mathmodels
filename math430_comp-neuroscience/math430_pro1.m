close all;
clear all;
% Hodgkin-Huxley model 

% Area is 1 cm^2 

I0 = 1.92:0.001:1.93; %initial current[nA]
A = 0; %SS Amplitude of wave [nA]
f = 9*10^-3; %SS frequency of current wave [Hz]
% D=0.02; %standard deviation for noise 
D = 0;
% tau=3; %[ms] time constant 

% m0 = 0.1; %initial Na activation m 
h0 = 0.1; %initial Na inactivation h
n0 = 0.5; %initial K activation n
% M0 = 0.9;

Vth = -20; %SS [mV] voltage threshold to be passed to be considered a spike 

Ena = 55; % GA sodium voltage [mV]55
Ek = -90; %GA potatssium [mV]
El = -80; %SS leak [mV]

cm = 1; %GA specific membrane capacitance[uF/cm^2] 

Gna = 24; % GA membrane conductance Na [mS/cm^2]
Gnap = 0.20; %SS persistent sodium current [mS/cm^2]
Gkdr = 3; %GA delayed rectifier potassium current 
Gm = 1; %SS muscarine conductance 

%calculate Gl -> membrane conductance 
Gl = 0.02; %GA leak current 

%time
T = 2200; %max time [msec]
dt = 0.03; % GA make sure dt is small enough 
t = 0:dt:T; 

ti = 100; %start of pulse [msec]
tf = 2100; %end of pulse [msec]
H = zeros(1, length(t)); %heaviside function 
H(floor(ti/dt):floor(tf/dt)) = 1;

% square pulse current 
Iapp = zeros(1, length(t)); %Current applied 

%initialize variables
V = zeros(1, length(t));
% m = zeros(1, length(t));
h = zeros(1, length(t));
n = zeros(1, length(t));
% M = zeros(1, length(t));
V(1) = El; %initial voltage 
% m(1) = m0;
h(1) = h0;
n(1) = n0;
% M(1) = M0;
spike = 0; %number of times the model spikes 
freq = zeros(1, length(I0));

% voltage dependent activation/inactiv. curves 
minf = @(v) 1/(1+exp(-(v+30)/9.5)); %GA 9.5 13.3
hinf = @(v) 1/(1+exp(-(v+53)/-7)); %GA 
% hinf = @(v) 1/(1+exp((v+82)/7)); %SS
% ninf = @(v) 1/(1+exp(-(v+30)/10)); %GA 
ninf = @(v) 1/(1+exp(-(v+35)/10)); %SS 
pinf = @(v) 1/(1 + exp(-(v+40)/5)); %GA persistent 

%voltage independent tau 
Tadj = 3.0^((36 - 22)/10); %SS 
Jh = @(v) (0.37 + 2.78/(1+exp(-(v+40.5)/-6))); %GA
% Jn = @(v) (0.37 + 1.85/(1+exp(-(v+27)/-15))); %GA 
Jn = @(v) (1000/3.3)/(Tadj * (exp((v+35)/40) + exp(-(v+35)/20)) ); %SS

%loop for each I0
for i=1:length(I0)
    for j=1:length(t)-1
        eta = randn;
%         Iapp(j+1) = (I0(i)+A*sin(2*pi*f*t(j))+sqrt(2*D*dt)*eta)* H(j); %amplitude of I0
        Iapp(j+1) = (I0(i)+A*sin(2*pi*f*t(j))+sqrt(2*D*dt)*eta)* H(j); %amplitude of I0
        
        %runge kutta 
        kv1 = (Iapp(j) - Gl*(V(j)-El) - Gna*((minf(V(j)))^3)*h(j)*(V(j)-Ena) - Gkdr*((n(j))^4)*(V(j)-Ek) - Gnap*pinf(V(j))*(V(j)-Ena) - Gm*n(j)*(V(j)-Ek))/cm; 
%         km1 = ( minf(V(j))-m(j) )/Jm(V(j));
        kh1 = ( hinf(V(j))-h(j) )/Jh(V(j));
        kn1 = ( ninf(V(j))-n(j) )/Jn(V(j));
%         kM1 = ( Minf(V(j))-M(j) )/JM(V(j));

        av1 = V(j) + kv1 * dt;
        am1 = minf(V(j));
        ah1 = h(j) + kh1 * dt;
        an1 = n(j) + kn1 * dt;
%         aM1 = M(j) + kn1 * dt;

        kv2 = (Iapp(j+1) - Gl*(av1-El) - Gna*(am1^3)*ah1*(av1-Ena) - Gkdr*(an1^4)*(av1-Ek) - Gnap*pinf(V(j))*(av1-Ena) - Gm*an1*(av1-Ek))/cm;
%         km2 = ( minf(av1)-am1 )/Jm(av1);
        kh2 = ( hinf(av1)-ah1 )/Jh(av1);
        kn2 = ( ninf(av1)-an1 )/Jn(av1);
%         kM2 = ( Minf(av1)-aM1 )/JM(av1);

        V(j+1) = V(j) + (kv1 + kv2)*dt/2; 
%         m(j+1) = m(j) + (km1 + km2)*dt/2;
        h(j+1) = h(j) + (kh1 + kh2)*dt/2;
        n(j+1) = n(j) + (kn1 + kn2)*dt/2;
%         M(j+1) = M(j) + (kM1 + kM2)*dt/2;
        
        %check if spike occurred 
        if V(j+1)>=Vth && t(j+1)>=(ti+500) %ignore the first 500ms 
            spike = spike + 1;
        end
    end
    
    disp(I0(i));
    disp(spike);
    
    %number of total spikes over duration of applied current 
    freq(i) = spike/(tf-ti);
    spike = 0;

    
    %only plot V vs t for I0=0,1,2,3,4 (whole numbers)
    if mod(I0(i),0.005)==0
%         hold on
        figure
        plot(t, V);
        axis([(ti+500) T -100 50]);
        xlabel('t (msec)');
        ylabel('V (mV)');
        title(['V vs t using Hudgkin Huxley Model where I0=',num2str(I0(i))])
       
        figure
        hold on
%         plot(t,1.95,'--');
        plot(t, Iapp);
        xlabel('t (msec)');
        ylabel('stimulus (nA)');
        axis([(ti+500) T 0 2.5]);        
        if I0(i) == 1.8
            title('Fluctuation-Driven');
        else 
            title('Mean-Driven');
        end

    end
end

figure 
plot(I0, freq);
xlabel('I0 (mV)');
ylabel('frequency (spike/msec)');
title('Freq vs I0 using Hodgkin Huxley Model');



