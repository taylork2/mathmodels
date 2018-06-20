%2017-09-18

% Passive membrane equation

clearvars;
close all

% Biophysical parameters

C = 1; %uF
El = -60;%mV
Ek = -85;%mV
Gl = (0.1);%S
Iapp = 1.8;%nA
D = 0;

Vth = -50; %mV
Vrst = -65; %mV
dG_sra = 0.1;
tau_sra = 10; %msec

% Time definitions

Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;

% Square wave (Heaviside function)

ti = dt;
tf = 10000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

% Initial conditions (for D = 0)

V = zeros(1,length(t));
V(1) = El;
W = zeros(1,length(t));
W(1) = 0;
% Computation of the solution (for D = 0)

for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)-W(j)*(V(j)-Ek)+Iapp*H(j))/C;
    kw1 = -W(j)/J_sra;
    av = V(j)+kv1*dt;
    aw = W(j)+kw1*dt;
    kv2 = (-Gl*(av-El)-aw*(av-Ek)+Iapp*H(j+1))/C;
    kw2 = -aw/J_sra;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;
    W(j+1) = W(j) + (kw1+kw2)*dt/2;
    
    if V(j+1) >= Vth
        V(j) = 60;
        V(j+1) = Vrst;
        
        W(j+1) = W(j)+ dG;
    end
end

tspk = zeros(1);
ISI = zeros(1);

% Initial conditions (for D > 0)

Vn = zeros(1,length(t));
Vn(1) = V(1);

% Computation of the solution (for D > 0)

for j=1:length(t)-1
    eta = randn;
    Vnaux = Vn(j)+sqrt(2*D*dt)*eta;
    kv1 = (-Gl*(Vnaux-El)+Iapp*H(j))/C;
    av = Vn(j)+kv1*dt;
    av = av + sqrt(2*D*dt)*eta;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    Vn(j+1) = Vn(j) + (kv1+kv2)*dt/2;
    Vn(j+1) = Vn(j+1) + sqrt(2*D*dt)*eta;   
end

% Graph of V(t)
figure
hold on
subplot(2,1,1)
plot(t,V,'b','linewidth',2);
xlabel('t')
ylabel('V');
set(gca,'fontsize',10);
axis([0 Tmax -80 80]);
%Graph of G_sra(t)
subplot(2,1,2)
plot(t,W,'r','linewidth',2);
xlabel('t')
ylabel('g_{sra}');
axis([0 Tmax  0 0.1]);
set(gca,'fontsize',10);
if D>0
    plot(t,Vn,'r','linewidth',2);
end

title('Voltage and G_{sra}')


