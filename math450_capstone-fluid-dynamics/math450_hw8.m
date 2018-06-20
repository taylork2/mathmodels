close all;
clearvars;
 
%initial variables 
N=50; %x steps 
m=10000; %number of time steps 
M=100000; %max iterations 
h0=1; %h initial
E=0.001; %episilon - amplitude of wave 
max = floor(N/4); %location of max of fxn 
 
p=1; %density
g=9.8; %gravity 
u=300; 
k=2*pi; 
o=0.0001; %surface tension
beta=p*g/(3*u); 
delta=o/(3*u); 
 
dx=1/N;
dt=1/m;
alpha = dt/(8*dx^2); 
gamma = dt/(8*dx^4);
 
x = linspace(0,1,N+1); %initialize x variables
h = h0 + E*sin(k*x); %initial values of h
hn=h; %create hn 
 
%to plot the slope 
y=zeros(1, m);
t=linspace(0,m,m);
 
plot(x, hn);
hold on
 
for j=1:M
    %numerically calculate next time step
    for i=2:N
        if i==2
            ghost1 = h(N);
            ghost2 = h(i+2);
        elseif i==N
            ghost1 = h(i-2);
            ghost2 = h(2);
        else 
            ghost1 = h(i-2);
            ghost2 = h(i+2);
        end
        hn(i) = h(i) + beta*alpha*( (h(i)+h(i+1))^3 * (h(i+1)-h(i)) - (h(i)+h(i-1))^3 * (h(i)-h(i-1))) + ...
           gamma*delta*( (h(i-1)+h(i))^3 * (ghost2-3*h(i+1)+3*h(i)-h(i-1)) - (h(i)+h(i-1))^3 * (h(i+1)-3*h(i)+3*h(i-1)-ghost1) );
    end
    
    %prevent each step from being graphed (for computational time's sake)
    if mod(j,1000)==0
        plot(x,hn);
        hold on
    end
    
    h=hn;
    
end

%graphing the wave 
figure(1)
title('Explicit solution for Height of water');
xlabel('x');
ylabel('h');
