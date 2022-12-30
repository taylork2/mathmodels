%2014-01-23

% Passive membrane equation

clear all
close all

C = 1;
Rm = 3333;
Gl = 0.3;
El = -50.6;
ENa = 55;
Ek = -75;
GNa = 120;
Gk = 36;

r = 0.0001*1; %Multiply radius by 10 times, to see spiking in 5th compartment
L = 0.159;
ra = 100;
rm = 3333;
c = 1;
%Iapp = 0.8:0.01:1;
Iapp = [0 1 2 3 4];
%El = -65;

A = 2*pi*r*L;
tau = c*rm
lambda = sqrt((r*rm)/(2*ra))

x = 0:1:1000;
dx = 0.001;

S1 = A*c;
g = (r*r^2)/(ra*L*(r^2*L +r^2*L))

%Functions



minf = @(V) 1/(1 + exp(-(V + 40)/9));
hinf = @(V) 1/(1 + exp((V + 62)/10));
ninf = @(V) 1/(1 + exp(-(V + 53)/16));
taum = @(V) 0.3;
tauh = @(V) 1 + 11/(1 + exp(V + 62)/10);
taun = @(V) 1 + 6/(1 + exp(V + 53)/16);

%Iapp = [0 1 2 3 4];


T = 4000;
dt = 0.1;
t = 0:dt:T;

ti = 1000;
tf = 3000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;
Vn1 = zeros(length(Iapp),length(t));
Vn2 = zeros(length(Iapp),length(t));
Vn3 = zeros(length(Iapp),length(t));
Vn4 = zeros(length(Iapp),length(t));
Vn5 = zeros(length(Iapp),length(t));
mf = zeros(length(Iapp),length(t));
hf = zeros(length(Iapp),length(t));
nf = zeros(length(Iapp),length(t));
V1 = zeros(1,length(t));
V1(1) = El;
V2 = zeros(1,length(t));
V2(1) = El;
V3 = zeros(1,length(t));
V3(1) = El;
V4 = zeros(1,length(t));
V4(1) = El;
V5 = zeros(1,length(t));
V5(1) = El;
m1 = zeros(1,length(t));
m2 = zeros(1,length(t));
m3 = zeros(1,length(t));
m4 = zeros(1,length(t));
m5 = zeros(1,length(t));
h1 = zeros(1,length(t));
h2 = zeros(1,length(t));
h3 = zeros(1,length(t));
h4 = zeros(1,length(t));
h5 = zeros(1,length(t));
n1 = zeros(1,length(t));
n2 = zeros(1,length(t));
n3 = zeros(1,length(t));
n4 = zeros(1,length(t));
n5 = zeros(1,length(t));
Vn1(:,1) = El;
Vn2(:,1) = El;
Vn3(:,1) = El;
Vn4(:,1) = El;
Vn5(:,1) = El;
m1(1) = 0;
m2(1) = 0;
m3(1) = 0;
m4(1) = 0;
m5(1) = 0;
h1(1) = 0;
h2(1) = 0;
h3(1) = 0;
h4(1) = 0;
h5(1) = 0;
n1(1) = 0;
n2(1) = 0;
n3(1) = 0;
n4(1) = 0;
n5(1) = 0;



spkcnt = 0;

for i=1:length(Iapp)
    for j=1:length(t)-1
        %Compartment 1 with injected current
        kv1 = (Iapp(i)*H(j) - Gl*(V1(j) - El) - GNa*m1(j)^3*h1(j)*(V1(j) - ENa) - Gk*n1(j)^4*(V1(j) - Ek))/C;
        km1 = (minf(V1(j)) - m1(j))/taum(V1(j));
        kh1 = (hinf(V1(j)) - h1(j))/tauh(V1(j));
        kn1 = (ninf(V1(j)) - n1(j))/taun(V1(j));
        av = V1(j) + kv1*dt;
        am = m1(j) + km1*dt;
        ah = h1(j) + kh1*dt;
        an = n1(j) + kn1*dt;
        kv2 = (Iapp(i)*H(j) - Gl*(av - El) - GNa*am^3*ah*(av - ENa) - Gk*an^4*(av - Ek))/C;
        km2 = (minf(av) - am)/taum(av);
        kh2 = (hinf(av) - ah)/tauh(av);
        kn2 = (ninf(av) - an)/taun(av);
        V1(j+1) = V1(j) + (kv1+kv2)*dt/2;
        m1(j+1) = m1(j) + (km1+km2)*dt/2;
        h1(j+1) = h1(j) + (kh1+kh2)*dt/2;
        n1(j+1) = n1(j) + (kn1+kn2)*dt/2;
        
        %Compartment 2 with NO injected current
        V2(j) = (V1(j)*g*rm)/(1 + rm*g);
                
        %Compartment 4 with NO injected current
        V4(j) = (V3(j)*g*rm)/(1 + rm*g);
        
        %Compartment 5 with NO injected current
        V5(j) = (V4(j)*g*rm)/(1 + rm*g);
        
        kv1 = (- Gl*(V5(j) - El) - GNa*m5(j)^3*h5(j)*(V5(j) - ENa) - Gk*n5(j)^4*(V5(j) - Ek))/C;
        km1 = (minf(V5(j)) - m5(j))/taum(V5(j));
        kh1 = (hinf(V5(j)) - h5(j))/tauh(V5(j));
        kn1 = (ninf(V5(j)) - n5(j))/taun(V5(j));
        av = V5(j) + kv1*dt;
        am = m5(j) + km1*dt;
        ah = h5(j) + kh1*dt;
        an = n5(j) + kn1*dt;
        kv2 = (- Gl*(av - El) - GNa*am^3*ah*(av - ENa) - Gk*an^4*(av - Ek))/C;
        km2 = (minf(av) - am)/taum(av);
        kh2 = (hinf(av) - ah)/tauh(av);
        kn2 = (ninf(av) - an)/taun(av);
        V5(j+1) = V5(j) + (kv1+kv2)*dt/2;
        m5(j+1) = m5(j) + (km1+km2)*dt/2;
        h5(j+1) = h5(j) + (kh1+kh2)*dt/2;
        n5(j+1) = n5(j) + (kn1+kn2)*dt/2;
        
    end
    Vn1(i,:) = V1;
    Vn2(i,:) = V2;
    Vn3(i,:) = V3;
    Vn4(i,:) = V4;
    Vn5(i,:) = V5;
    %mf(i,:) = m;
    %hf(i,:) = h;
    %nf(i,:) = n;
end



freq = zeros(length(Iapp));


for i = 1:length(freq)
    spikethreshold = 0;
    [spikeamp, idx] = findpeaks(Vn5(i,:)); 
    xx = find(spikeamp > spikethreshold); 
        spikeamp = spikeamp(xx); 
        idx = idx(xx);
    ifreq = diff(t(idx));
    ifreq = 1000 ./ifreq;
    if(isempty(ifreq))
        ifreq = zeros(1,24);
    end

    freq(i) = ifreq(1);
end

freq = freq(:,1)




figure

hold on


subplot(5,1,1)
plot(t,Vn1(5,:),'r','linewidth',2);
subplot(5,1,2)
plot(t,Vn2(5,:),'r','linewidth',2);
subplot(5,1,3)
plot(t,Vn3(5,:),'r','linewidth',2);
subplot(5,1,4)
plot(t,Vn4(5,:),'r','linewidth',2);
subplot(5,1,5)
plot(t,Vn5(5,:),'r','linewidth',2);



%plot(t,Vn(3,:),'r','linewidth',2);
%plot(Iapp,freq,'*')
%axis([0 10000 -75 -60]);

