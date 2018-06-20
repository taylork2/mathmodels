close all;
clear all;
% Hodgkin-Huxley model 

% Area is 1 cm^2 

I0 = 0:.25:4; %initial current[mV]
m0 = 0.1; %initial Na activation m 
h0 = 0.9; %initial Na inactivation h
n0 = 0.75; %initial K activation n

Vth = -40; %voltage threshold to be passed to be considered a spike 

cm = 1; %specific membrane capacitance[uF/cm^2] 
rm = 3333*10^-3; %specific membrane resistance[mohm*cm^2] 
El = -50.6; %potential [mV] 
Ena = 55; %potential Na [mV]
Ek = -75; %potential K [mV]
Gna = 120; % membrane conductance Na [mS/cm^2]
Gk = 36; % membrane conductance K [mS/cm^2]

%calculate Gl -> membrane conductance 
Gl = 1/rm;

%time
T = 300; %max time [msec]
dt = 0.1; %make sure dt is small enough 
t = 0:dt:T; 

ti = 200; %start of pulse [msec]
tf = 2200; %end of pulse [msec]
H = zeros(1, length(t)); %heaviside function 
H(floor(ti/dt):floor(tf/dt)) = 1;

% square pulse current 
Iapp = zeros(1, length(t)); %Current applied 

%initialize variables
V = zeros(1, length(t));
m = zeros(1, length(t));
h = zeros(1, length(t));
n = zeros(1, length(t));
V(1) = El; %initial voltage 
m(1) = m0;
h(1) = h0;
n(1) = n0;
spike = 0; %number of times the model spikes 
freq = zeros(1, length(I0));

% voltage dependent activation/inactiv. curves 
minf = @(v) 1/(1+exp(-(v+40)/9));
hinf = @(v) 1/(1+exp((v+62)/10));
ninf = @(v) 1/(1+exp(-(v+53)/16));

%voltage independent tau 
Jm = @(v) 0.3;
Jh = @(v) (12 + exp((v+62)/10))/(1 + exp((v+62)/10));
Jn = @(v) (7 + exp((v+53)/16))/(1 + exp((v+53)/16));

%loop for each I0
for i=1:length(I0)
    for j=1:length(t)-1
        Iapp(j+1) = I0(i) * H(j); %amplitude of I0
        
        %runge kutta 
        kv1 = (Iapp(j) - Gl*(V(j)-El) - Gna*((m(j))^3)*h(j)*(V(j)-Ena) - Gk*((n(j))^4)*(V(j)-Ek))/cm; 
        km1 = ( minf(V(j))-m(j) )/Jm(V(j));
        kh1 = ( hinf(V(j))-h(j) )/Jh(V(j));
        kn1 = ( ninf(V(j))-n(j) )/Jn(V(j));

        av1 = V(j) + kv1 * dt;
        am1 = m(j) + km1 * dt;
        ah1 = h(j) + kh1 * dt;
        an1 = n(j) + kn1 * dt;

        kv2 = (Iapp(j+1) - Gl*(av1-El) - Gna*(am1^3)*ah1*(av1-Ena) - Gk*(an1^4)*(av1-Ek))/cm;
        km2 = ( minf(av1)-am1 )/Jm(av1);
        kh2 = ( hinf(av1)-ah1 )/Jh(av1);
        kn2 = ( ninf(av1)-an1 )/Jn(av1);

        V(j+1) = V(j) + (kv1 + kv2)*dt/2; 
        m(j+1) = m(j) + (km1 + km2)*dt/2;
        h(j+1) = h(j) + (kh1 + kh2)*dt/2;
        n(j+1) = n(j) + (kn1 + kn2)*dt/2;
        
        %check if spike occurred 
        if V(j+1)>=Vth && t(j+1)>=ti
            spike = spike + 1;
        end
    end
    
    %number of total spikes over duration of applied current 
    freq(i) = spike/(tf-ti);
    
    %only plot V vs t for I0=0,1,2,3,4 (whole numbers)
    if mod(I0(i),1)==0
        figure
        plot(t, V);
        axis([40 T -100 40]);
        xlabel('t (msec)');
        ylabel('V (mV)');
        title(['V vs t using Hudgkin Huxley Model where I0=',num2str(I0(i))])
    end
end

figure 
plot(I0, freq);
xlabel('I0 (mV)');
ylabel('frequency (spike/msec)');
title('Freq vs I0 using Hodgkin Huxley Model');



