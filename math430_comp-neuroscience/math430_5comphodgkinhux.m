close all; 

I0 = 0.008; %initial current[mA]
m0 = 0.1; %initial Na activation m 
h0 = 0.9; %initial Na inactivation h
n0 = 0.75; %initial K activation n

Vth = -40; %voltage threshold to be passed to be considered a spike 

Ena = 55; %potential Na [mV]
Ek = -75; %potential K [mV]
Gna = 120; % membrane conductance Na [mS/cm^2]
Gk = 36; % membrane conductance K [mS/cm^2]

cm = 1; %capacitance [uF/cm^2]
rm = 3333*10^-3; % membrane resistance [ohm->mohm cm^2]
rl = 100*10^-3; % cytoplasm resistance [ohm->mohm cm]
El = -65; % [mV->V]
p = 1*10^-4; %radius [um->cm]
p = 1*10^-3;
L0 = 0.159; %length of compartment [cm]
L0 = 0.09;

N = 5; %num of compartments 
L = [L0 L0 L0 L0 L0]; %each compartment same length
a = [p p p p p]; %each comp same radius
c = [cm cm cm cm cm]; %each comp same capacitance 
R = [rm rm rm rm rm]; %each comp same resistance
V0 = El;

%calculate Gl -> membrane conductance 
Gl = 1/rm;

%time 
T = 100; %max time [s]
dt = 0.01;
t = 1:dt:T;

ti = 20; %start of current [s]
tf = 80; %end of current [s]
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

% Initialize m,n and h 
m = zeros(N, length(t));
m_k1 = zeros(1,N);
m_a = zeros(1,N);
m_k2 = zeros(1,N);
h = zeros(N, length(t));
h_k1 = zeros(1,N);
h_a = zeros(1,N);
h_k2 = zeros(1,N);
n = zeros(N, length(t));
n_k1 = zeros(1,N);
n_a = zeros(1,N);
n_k2 = zeros(1,N);

%initialize V for each compartment 
V = zeros(N, length(t));
V_k1 = zeros(1, N);
V_av = zeros(1, N);
V_k2 = zeros(1, N);

%set V0 for each V array 
for i=1:N
    V(i,1) = V0;
    m(i,1) = m0;
    h(i,1) = h0;
    n(i,1) = n0;
end

% voltage dependent activation/inactiv. curves 
minf = @(v) 1/(1+exp(-(v+40)/9));
hinf = @(v) 1/(1+exp((v+62)/10));
ninf = @(v) 1/(1+exp(-(v+53)/16));

%voltage independent tau 
Jm = @(v) 0.3;
Jh = @(v) (12 + exp((v+62)/10))/(1 + exp((v+62)/10));
Jn = @(v) (7 + exp((v+53)/16))/(1 + exp((v+53)/16));

g = zeros(N, N); %g(1,2) coupling resistance b/w comp 1 and 2
A = zeros(1,N); %Area 

%calc g [ohm^-1 cm^-2]
for i=1:N
    for j=1:N
        g(i,j) = (a(i)*a(j)^2) / (rl*L(i)*(a(j)^2*L(i) + a(i)^2*L(j)));
    end
end

%calc Areas of comps [cm^2]
for i=1:N
    A(i) = L(i)*2*pi*a(i);
end

Iapp = zeros(N, length(t)); %Current applied 

for k=1:length(t)-1
    % Only apply current to first compartment
    Iapp(1,k+1) = I0 * H(k); %amplitude of I0

    % calculate the V_k1 
    for i=1:N
        if i==1
            Vb = 0;
            Va = (V(2,k)-V(1,k)).*g(1,2); 
        elseif i==N
            Vb = (V(N-1,k)-V(N,k)).*g(N,N-1);
            Va = 0;
        else
            Vb = (V(i-1,k)-V(i,k)).*g(i,i-1);
            Va = (V(i+1,k)-V(i,k)).*g(i+1,i);
        end
        V_k1(i) = 1/cm*((Vb + Va) + Iapp(i,k)/A(i) + (- Gl*(V(i,k)-El) - Gna*((m(i,k))^3)*h(i,k)*(V(i,k)-Ena) - Gk*((n(i,k))^4)*(V(i,k)-Ek))/rm); 
        m_k1(i) = ( minf(V(i,k))-m(i,k) )/Jm(V(i,k));
        h_k1(i) = ( hinf(V(i,k))-h(i,k) )/Jh(V(i,k));
        n_k1(i) = ( ninf(V(i,k))-n(i,k) )/Jn(V(i,k));
    end

    % Calculate V_av 
    for i=1:N
        V_av(i) = V(i,k) + V_k1(i)* dt;
        m_a(i) = m(i,k) + m_k1(i) * dt;
        h_a(i) = h(i,k) + h_k1(i) * dt;
        n_a(i) = n(i,k) + n_k1(i) * dt;
    end
    
    % Calculate V_k2 
    for i=1:N
        if i==1
            Vb = 0;
            Va = (V_av(2)-V_av(1)).*g(1,2); 
        elseif i==N
            Vb = (V_av(N-1)-V_av(N)).*g(N,N-1);
            Va = 0;
        else
            Vb = (V_av(i-1)-V_av(i)).*g(i,i-1);
            Va = (V_av(i+1)-V_av(i)).*g(i+1,i);
        end
        %(Vb + Va) + Iapp(i,k)/A(i) - 
        V_k2(i) = 1/cm*((Vb + Va) + Iapp(i,k)/A(i) + (- Gl*(V_av(i)-El) - Gna*((m_a(i))^3)*h_a(i)*(V_av(i)-Ena) - Gk*((n_a(i))^4)*(V_av(i)-Ek))/rm); 
        m_k2(i) = ( minf(V_av(i))-m_a(i) )/Jm(V_av(i));
        h_k2(i) = ( hinf(V_av(i))-h_a(i) )/Jh(V_av(i));
        n_k2(i) = ( ninf(V_av(i))-n_a(i) )/Jn(V_av(i));
    end
    
    %calc V 
    for i=1:N
        V(i,k+1) = V(i,k) + (V_k1(i) + V_k2(i))*dt/2;
        m(i,k+1) = m(i,k) + (m_k1(i) + m_k2(i))*dt/2;
        h(i,k+1) = h(i,k) + (h_k1(i) + h_k2(i))*dt/2;
        n(i,k+1) = n(i,k) + (n_k1(i) + n_k2(i))*dt/2;
    end
end

figure 
hold on
subplot(4,1,4);
plot(t, Iapp(1,:));

subplot(4,1,[1 2 3]);
hold on
plot(t, V(1,:), 'r', 'linewidth', 2);   
plot(t, V(2,:));
plot(t, V(3,:));
plot(t, V(4,:));
plot(t, V(5,:), 'linewidth', 2);
legend("Comp 1", "Comp 2", "Comp 3", "Comp 4", "Comp 5");
xlabel('t(ms)');
ylabel('V(mV)');
title('V vs t with change in length & radius');


