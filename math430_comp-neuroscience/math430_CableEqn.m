%Cable equation/compartmental model
clear all;
close all;

% 1A)
%% 
p = 2.5*10^-4; % radius[um->cm]
rl = 200*10^-3; % Resistance of cytoplasm [ohm->kohm cm]
rm = 20000*10^-3; % specific membrane resistance [ohm->kohm cm^2]
Cm = 1; % specific membrane capacitance [uF/cm^2]

Jm = Cm*rm; 

lamda = sqrt(p*rm/(2*rl));

disp(['Time constant Jm = ', num2str(Jm)]);
disp(['Length constant lamda = ', num2str(lamda)]);


%% 
%1B) 
clear all;
close all;
d = 4*10^-4; %diameter [um->cm]
a = d/2; % radius [cm]
rl = 100*10^-3; % Resistance of cytoplasm [ohm->kohm cm]
rm = 10000*10^-3; % specific membrane resistance [ohm->kohm cm^2]
Iapp = 0.25; % Applied current [nA]
Vrest = -65; % voltage [mV]

X = 10; %max length [cm]
dx = 0.01; 
x = 0:dx:X;

%calculate lamda length constant 
lam = sqrt(a*rm/(2*rl)); %[cm]

%steady state equation 
Vs = (lam*rl/(pi*a^2)*Iapp*exp(-x/lam)); %[mV]

V = Vs - Vrest; 
V0 = V(1); %max value 
for i=1:length(V)
    if (V(i)/V0 <= 0.63)
        disp(i*dx);
        break;
    elseif (i==length(V))
        disp("There are no values of x where function decreases by 63% of its max value");
    end
end   

figure(1)
plot(x, Vs);
xlabel('x (cm)');
ylabel('V (mV)');
title('Cable Equation V_{s} vs x');

figure(2)
plot(x, V);
xlabel('x (cm)');
ylabel('V (mV)');
title('Cable Equation V_{s}-V_{rest} vs x');


% 2)
%% 
clear all;
close all; 

cm = 1*10-4; %capacitance [uF->F/cm^2]
rm = 3333*10^-3; % membrane resistance [ohm cm^2]
rl = 100*10^-3; % cytoplasm resistance [ohm cm]
El = -65; % [mV->V]

% 2A)
Jm = cm*rm; %[sec]
disp(['Time constant Jm = ', num2str(Jm)]);

% 2B)
p = 1*10^-4; %radius [um->cm]
L0 = 0.159; %length of compartment [cm]
A0 = 2*pi*p*L0; %[cm^2]
disp(['The product A*Cm is = ', num2str(A0*cm)]); %[F cm]

% 2C) THIS IS g!!!
Djk = (p*p^2) / (rl*L0*(p^2*L0 + p^2*L0));
disp(['D_{jk}/C_{m} = ', num2str(Djk/cm)]);

% 2D) 
% N = 5; %num of compartments 
% L = [L0 L0 L0 L0 L0]; %each compartment same length
% a = [p p p p p]; %each comp same radius
% c = [cm cm cm cm cm]; %each comp same capacitance 
% R = [rm rm rm rm rm]; %each comp same resistance
% V0 = El;
% 
% gl = 1/rl; %S
% 
% %time 
% T = 10000; %max time [s]
% dt = 0.01;
% t = 1:dt:T;
% 
% ti = 100; %start of current [s]
% tf = 7000; %end of current [s]
% H = zeros(1,length(t));
% H(floor(ti/dt):floor(tf/dt))=1;
% 
% g = zeros(N-1, N-1); %g(1,2) coupling resistance b/w comp 1 and 2
% r = zeros(N-1, N-1);
% A = zeros(1,N); %Area 
% 
% d = 2*p; %diameter [cm]
% lamda = sqrt(d*rm/(4*rl));
% 
% %calc g [ohm^-1 cm^-2]
% for i=1:N
%     for j=1:N
%         g(i,j) = (a(i)*a(j)^2) / (rl*L(i)*(a(j)^2*L(i) + a(i)^2*L(j)));
%     end
% end
% 
% %calc Areas of comps [cm^2]
% for i=1:N
%     A(i) = L(i)*2*pi*a(i);
% end
% 
% %initialize V for each compartment 
% V = zeros(N, length(t));
% V_k1 = zeros(N, length(t));
% V_av = zeros(1, N);
% V_k2 = zeros(N, length(t));
% 
% %set V0 for each V array 
% for i=1:N
%     V(i,1) = V0;
% end
% 
% Iel0 = 0; %[A]
% Iel_app = 1;
% 
% ie = zeros(N, length(t));
% Iel = Iel0 * ones(1, N); %electrode current of each comp
% 
% %only inject current into First compartment
% for j=1:length(t)
%     ie(1,j) = H(j)*Iel_app/A(1);
% end
% 
% 
% for k=1:length(t)-1
%     
%     % calculate the V_k1 
%     % Use passive membrane equation for initial compartment
%     V_k1(1,k) = (1/c(1))*(-gl*(V(1,k)-El)+ie(1,k)*H(k));
%     
%     % Calculate remaining compartments using compartmental model 
%     V_k1(2,k) = (1/c(2)).*(( (V(1,k)-V(2,k)).*g(2,1)) + (V(3,k)-V(2,k)).*g(2,3) - (V(2,k)./R(2))+ie(2,k));
%     V_k1(3,k) = (1/c(3)).*(( (V(2,k)-V(3,k)).*g(3,2)) - (V(3,k)./R(3)) + ie(3,k) );
%     V_k1(4,k) = (1/c(4)).*(( (V(1,k)-V(4,k)).*g(4,1)) + (V(5,k)-V(4,k)).*g(4,5) - (V(4,k)./R(4))+ie(4,k) );
%     V_k1(5,k) = (1/c(5)).*(( (V(5,k)-V(5,k)).*g(4,5)) - (V(5,k)./R(5)) + ie(5,k) );
% 
%     % Calculate V_av 
%     for i=1:N
%         V_av(i) = V(i,k) + V_k1(i,k)*dt;
%     end
%     
%     %calc V_k2 [V]
%     % Passive membrane equation 
%     V_k2(1,k) = (1/c(1))*(-gl*(V_av(1)-El)+ie(1,k+1)*H(k+1));
%     
%     V_k2(2,k) = (1/c(2)).*(( (V_av(1) - V_av(2) ).*g(2,1)) + (V_av(3) - V_av(2) ).*g(2,3) - (V_av(2)./R(2))+ie(2,k));
%     V_k2(3,k) = (1/c(3)).*(( (V_av(2) - V_av(3) ).*g(3,2)) - (V_av(3)/R(3)) + ie(1,k));
%     V_k2(4,k) = (1/c(4)).*(( (V_av(1) - V_av(4) ).*g(4,1)) + (V_av(5)-V_av(4)).*g(4,5) - (V_av(4)./R(4)) + ie(4,k));
%     V_k2(5,k) = (1/c(5)).*(( (V_av(4) - V_av(5) ).*g(5,4)) - (V_av(5)./R(5)) + ie(5,k));            
%     
%     
%     %calc V 
%     for i=1:N
%         V(i,k+1) = V(i,k) + (V_k1(i,k) + V_k2(i,k))*dt/2;
%     end
% end
% 
% figure 
% hold on
% %plot(t, ie(1,:));
% plot(t, V(1,:));   
% plot(t, V(2,:));
% plot(t, V(3,:));
% plot(t, V(4,:));
% plot(t, V(5,:));
% legend("1", "2", "3", "4", "5");
% xlabel('t(s)');
% ylabel('V(V)');
% title('V vs t');



    









