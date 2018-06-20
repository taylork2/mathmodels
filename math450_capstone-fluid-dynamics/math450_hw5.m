close all;
clearvars;

%initial variables 
N=50; %x steps 
m=10000; %number of time steps 
h0=1; %h initial
E=0.001; %episilon - amplitude of wave 
max = floor(N/4); %location of max of fxn 

p=1;
g=9.8; 
u=30;
k=2*pi;
beta=p*g/(3*u); %will = 1/9 for testing 

dx=1/N;
dt=1/m;
alpha = dt/(8*dx^2); 

x = linspace(0,1,N+1); %initialize x variables
h = h0 + E*cos(k*x); %initial values of h
hn=h; %create hn 

%to plot the slope 
y=zeros(1, m);
t=linspace(0,m,m);

%plot(x, hn);
%hold on

for j=1:m
    %numerically calculate next time step
    for i=2:N
        hn(i) = h(i) + beta*alpha*( (h(i)+h(i+1))^3 * (h(i+1)-h(i)) - (h(i)+h(i-1))^3 * (h(i)-h(i-1)));
    end
    
    %prevent each step from being graphed (for computational time's sake)
    if mod(j,50)==0
        plot(x,hn);
        hold on
    end
    
    %to get the slope 
    y(j)=log(abs((h(max)-h0)/(E*h0)));
    
    h=hn;
    
end

%slope 
w = -beta * h0^3 * k^2 / m
coef = polyfit(t, y, 1);
slope = coef(1)

%graphing the wave 
figure(1)
title('Numerical solution for Height of water');
xlabel('x');
ylabel('h');

%graphing the slope 
figure(2)
plot(t, y);
title('Slope at x_{max} vs t');
xlabel('t');
ylabel('ln|(h(x_{max},t)-h0)/E*h0)|');


