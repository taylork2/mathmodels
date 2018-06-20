close all;
clearvars;

%initial variables 
N=128; %x steps 
m=100000; %number of time steps 
T=10; %max Time iterations 
M=100; %max Newton's method iterations  
h0=1; %h initial 
E=0.01; %episilon - amplitude of wave 
max = floor(N/4); %location of max of fxn 

%to plot the slope 
y=zeros(1, T);
time=linspace(0,T/m,T);

p=1; %density
g=9.8; %gravity 
u=1; 
k=2*pi; 
o=0.01; %surface tension
beta=p*g/(3*u); 
delta=o/(3*u); 

dx=1/N;
dt=1/m;
alpha = dt/(8*dx^2); 
gamma = dt/(8*dx^4);

x = linspace(0,1,N+1); %initialize x variables
h_plot = (h0 + E*sin(k*x))'; %initial values of h
hor = h_plot;
h = h_plot(1:N);

% hn_plot=h_plot; %create hn 
hn=h;

oldh=h;

plot(x, h_plot);
hold on

error = zeros(N,1);
conv = ones(N,1) * 10^-7; %value to determine convergence  
maxval = 10^3; %value to determine divergence 

%Jacobian stuff 
s=alpha*beta;
v=gamma*delta;
% v=0;

% df1 = @(h1, h2, h3, h4, h5) -v*(h3+h2)^3; 
% 
% df2 = @(h1, h2, h3, h4, h5) s*( 3*(h3+h2)^2 * (h3-h2) - (h3+h2)^3 ) + ...
%     v*( 3*(h3+h2)^3 - 3*(h3+h2)^2*(h4-3*h3+3*h2-h1) + (h4+h3)^3); 
% 
% df3 = @(h1, h2, h3, h4, h5) 1 - s*( 3*( h4+h3 )^2*( h4-h3 ) - (h4+h3)^3 - ...
%     3*(h3+h2)^2*(h3-h2) - (h3+h2)^3 ) - v*( 3*(h4+h3)^3 + 3*(h4+h3)^2*(h5-3*h4+3*h3-h2) - ... 
%     ( -3*(h3+h2)^3 + 3*(h3+h2)^2*(h4-3*h3+3*h2-h1) ));
% 
% df4 = @(h1, h2, h3, h4, h5) -s*( 3*(h4+h3)^2 * (h4-h3) + (h4+h3)^3 ) - ...
%     v * ( -3*(h4+h3)^3 + 3*(h4+h3)^2*(h5-3*h4+3*h3-h2) - (h3+h2)^3); 
% 
% df5 = @(h1, h2, h3, h4, h5) -v*(h4+h3)^3; 

S = (dt)/(24*u*(dx^2));
a1 = p*g;       %constant infront of gravity scheme
a2 = o/(dx^2); %contant in front of surface tension scheme

%derivative wrt x1
df1=@(x1,x2,x3,x4,x5) -S*a2*(x2 + x3)^3;  
%derivative wrt x2
df2=@(x1,x2,x3,x4,x5) S*(a2*(3*(x2 + x3)^3 - 3*(x2 + x3)^2*(x1 - 3*x2 + 3*x3 - x4) + ...    
        (x3 + x4)^3) - a1*((x2 + x3)^3 + 3*(x2 + x3)^2*(x2 - x3)));   
%derivative wrt x3    
df3=@(x1,x2,x3,x4,x5) 1 - S*(a2*(3*(x2 + x3)^2*(x1 - 3*x2 + 3*x3 - x4) - ...      
        3*(x3 + x4)^2*(x2 - 3*x3 + 3*x4 - x5) + 3*(x2 + x3)^3 + 3*(x3 + x4)^3) - ...
        a1*((x2 + x3)^3 + (x3 + x4)^3 - 3*(x2 + x3)^2*(x2 - x3) + 3*(x3 + x4)^2*(x3 - x4)));  
%derivative wrt x4    
df4=@(x1,x2,x3,x4,x5) S*(a2*(3*(x3 + x4)^2*(x2 - 3*x3 + 3*x4 - x5) + ...
        (x2 + x3)^3 + 3*(x3 + x4)^3) - a1*((x3 + x4)^3 - 3*(x3 + x4)^2*(x3 - x4)));
%derivative wrt x5    
df5=@(x1,x2,x3,x4,x5) -S*a2*(x3 + x4)^3;

J=zeros(N,N);

%Gravity scheme (original)
g1 = @(x2,x3,x4) ((x4+x3)^3) *(x4-x3) - ((x3+x2)^3)*(x3-x2);
%Surface Tension scheme
g2 = @(x1,x2,x3,x4,x5) ((x4+x3)^3) * (x5-3*x4+3*x3-x2) - ((x3+x2)^3) * (x4-3*x3+3*x2-x1);
%Total scheme with Gravity and Surface Tension
f= @(oldx, x1,x2,x3,x4,x5) x3-(s*g1(x2,x3,x4) + v*g2(x1,x2,x3,x4,x5) ) - oldx;

% f = @(x, h1, h2, h3, h4, h5) h3 - x - s * ( (h4+h3)^3 * (h4-h3) - (h3+h2)^3*(h3-h2)) - ...
%     v * ((h4+h3)^3*(h5- 3*h4 + 3*h3 -h2) - (h3+h2)^3*(h4- 3*h3 + 3*h2 -h1));
b = h; %first guess

 
for t=1:T
    
    for i=1:M
        %h(N) and h(N-1) wrap around bc ring topology 
        b(1)=f(oldh(1),h(N-1),h(N),h(1),h(2),h(3));
        b(2)=f(oldh(2),h(N)  ,h(1),h(2),h(3),h(4));
        %h(1) 
        b(N)   = f(oldh(N),  h(N-2),h(N-1),h(N)  ,h(1),h(2));
        b(N-1) = f(oldh(N-1),h(N-3),h(N-2),h(N-1),h(N),h(1));

        for j=3:N-2
            b(j) = f(oldh(j),h(j-2),h(j-1),h(j),h(j+1),h(j+2));
        end
      
        %Create Jacobian 
        J(1,1) = df3(h(N-1),h(N), h(1), h(2), h(3)); 
        J(1,2) = df4(h(N-1),h(N), h(1), h(2), h(3));
        J(1,3) = df5(h(N-1),h(N), h(1), h(2), h(3));
        J(1,N-1) = df1(h(N-1),h(N), h(1), h(2), h(3));
        J(1,N) = df2(h(N-1),h(N), h(1), h(2), h(3));
        
        J(2,1) = df2(h(N), h(1), h(2), h(3), h(4));
        J(2,2) = df3(h(N), h(1), h(2), h(3), h(4));
        J(2,3) = df4(h(N), h(1), h(2), h(3), h(4));
        J(2,4) = df5(h(N), h(1), h(2), h(3), h(4));
        J(2,N) = df1(h(N), h(1), h(2), h(3), h(4));

        J(N-1,N-3) = df1(h(N-3),h(N-2), h(N-1), h(N), h(1));
        J(N-1,N-2) = df2(h(N-3),h(N-2), h(N-1), h(N), h(1));
        J(N-1,N-1) = df3(h(N-3),h(N-2), h(N-1), h(N), h(1));
        J(N-1,N)   = df4(h(N-3),h(N-2), h(N-1), h(N), h(1));
        J(N-1,1)   = df5(h(N-3),h(N-2), h(N-1), h(N), h(1));

        J(N,N-2) = df1(h(N-2), h(N-1), h(N), h(1), h(2));
        J(N,N-1) = df2(h(N-2), h(N-1), h(N), h(1), h(2));
        J(N,N)   = df3(h(N-2), h(N-1), h(N), h(1), h(2));
        J(N,1)   = df4(h(N-2), h(N-1), h(N), h(1), h(2));
        J(N,2)   = df5(h(N-2), h(N-1), h(N), h(1), h(2));
        
        for j=3:N-2
            J(j,j-2) = df1(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j-1) = df2(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j  ) = df3(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j+1) = df4(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j+2) = df5(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
        end
        
        hn = h - J'\b; 
        error = abs(hn-h);
        if error < conv
            h=hn;
            break;
        elseif abs(hn) > maxval
            disp("diverge");
            break;
        end     
        h=hn;
    end
    
    oldh = h;

    %prevent each step from being graphed (for computational time's sake)
%     if mod(t,50)==0
%     if t == T
        h_plot(1:N) = h;
        h_plot(N+1) = h(1);
        plot(x,h_plot);
%     end

    %to get the slope 
    y(t)=log(abs((h(max)-h0)/(E*h0)));
end

w = -(h0).^3*k^2*(o*k^2 + p*g)/(3*u*m)
% w = -beta * h0^3 * k^2 / m
coef = polyfit(time, y, 1);
slope = coef(1)

%graphing the wave 
figure(1)
title('Implicit Height of water with Surface Tension and periodic boundaries');
xlabel('x');
ylabel('height');

figure(2)
plot(time, y);
title('Slope at x_{max} vs t');
xlabel('t');
ylabel('ln|(h(x_{max},t)-h0)/E*h0)|');
