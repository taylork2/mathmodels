close all;
clearvars;

%initial variables 
N=128/2; %x steps (Used 128 for resolution for Basilisk) 
m=1/.003043; %number of time steps .003043
T=738; %max Time iterations (number of points from Basilisk)
M=100; %max Newton's method iterations  
h0=.1; %h initial 
E=0.01; %episilon - amplitude of wave 
max = floor(N/4); %location of max of fxn 

dx=1/N;
dt=1/m;
% dt=3.04290519077*10^-3;
C1 = dt/(8*dx^2); 
C2 = dt/(8*dx^4);

%to plot the slope 
y=zeros(1, T);
% basdt = 3.04290519077*10^-3; %dt used in Basilisk
time=linspace(0,T*dt,T); 
hmax=y;


p=1; %density
g=0; %gravity 
u=1; 
k=2*pi; 
o=1; %surface tension
c1=p*g/(3*u); 
c2=-o/(3*u); 
% c1 = 0; %no gravity 


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
S1=C1*c1;
S2=C2*c2;

df1 = @(h1, h2, h3, h4, h5) -S2*(h3+h2)^3; 

df2 = @(h1, h2, h3, h4, h5) S1*( 3*(h3+h2)^2 * (h3-h2) - (h3+h2)^3 ) + ...
    S2*( 3*(h3+h2)^3 + 3*(h3+h2)^2*(h4-3*h3+3*h2-h1) + (h4+h3)^3); 

df3 = @(h1, h2, h3, h4, h5) 1 - S1*( 3*( h4+h3 )^2*( h4-h3 ) - (h4+h3)^3 - ...
    3*(h3+h2)^2*(h3-h2) - (h3+h2)^3 ) - S2*( 3*(h4+h3)^3 + 3*(h4+h3)^2*(h5-3*h4+3*h3-h2) - ... 
    ( -3*(h3+h2)^3 + 3*(h3+h2)^2*(h4-3*h3+3*h2-h1) ));

df4 = @(h1, h2, h3, h4, h5) -S1*( 3*(h4+h3)^2 * (h4-h3) + (h4+h3)^3 ) - ...
    S2 * ( -3*(h4+h3)^3 + 3*(h4+h3)^2*(h5-3*h4+3*h3-h2) - (h3+h2)^3); 

df5 = @(h1, h2, h3, h4, h5) -S2*(h4+h3)^3; 

% %derivative wrt x1
% df1=@(x1,x2,x3,x4,x5) -S2*(x2 + x3)^3;  
% %derivative wrt x2
% df2=@(x1,x2,x3,x4,x5) (S2*(3*(x2 + x3)^3 - 3*(x2 + x3)^2*(x1 - 3*x2 + 3*x3 - x4) + ...    
%         (x3 + x4)^3) - S1*((x2 + x3)^3 + 3*(x2 + x3)^2*(x2 - x3)));   
% %derivative wrt x3    
% df3=@(x1,x2,x3,x4,x5) 1 - (S2*(3*(x2 + x3)^2*(x1 - 3*x2 + 3*x3 - x4) - ...      
%         3*(x3 + x4)^2*(x2 - 3*x3 + 3*x4 - x5) + 3*(x2 + x3)^3 + 3*(x3 + x4)^3) - ...
%         S1*((x2 + x3)^3 + (x3 + x4)^3 - 3*(x2 + x3)^2*(x2 - x3) + 3*(x3 + x4)^2*(x3 - x4)));  
% %derivative wrt x4    
% df4=@(x1,x2,x3,x4,x5) (S2*(3*(x3 + x4)^2*(x2 - 3*x3 + 3*x4 - x5) + ...
%         (x2 + x3)^3 + 3*(x3 + x4)^3) - S1*((x3 + x4)^3 - 3*(x3 + x4)^2*(x3 - x4)));
% %derivative wrt x5    
% df5=@(x1,x2,x3,x4,x5) -S2*(x3 + x4)^3;

J=zeros(N,N);

%Gravity scheme (original)
g1 = @(x2,x3,x4) ((x4+x3)^3) *(x4-x3) - ((x3+x2)^3)*(x3-x2);
%Surface Tension scheme
g2 = @(x1,x2,x3,x4,x5) ((x4+x3)^3) * (x5-3*x4+3*x3-x2) - ((x3+x2)^3) * (x4-3*x3+3*x2-x1);
%Total scheme with Gravity and Surface Tension
% f= @(oldx, x1,x2,x3,x4,x5) x3-(S1*g1(x2,x3,x4) + S2*g2(x1,x2,x3,x4,x5) ) - oldx;

f = @(x, h1, h2, h3, h4, h5) h3 - x - S1 * ( (h4+h3)^3 * (h4-h3) - (h3+h2)^3*(h3-h2)) - ...
    S2 * ((h4+h3)^3*(h5- 3*h4 + 3*h3 -h2) - (h3+h2)^3*(h4- 3*h3 + 3*h2 -h1));
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
    hmax(t) = h(max)-h0;
    y(t)=log(abs((h(max)-h0)/(E*h0)));
end

coef = polyfit(time, y, 1);
slope = coef(1)

%graphing the wave 
figure(1)
title('Implicit Height of water with Surface Tension and periodic boundaries');
xlabel('x');
ylabel('height');

figure(2)
plot(time, y);
title('Linear Stability at x_{max} vs t');
xlabel('t');
ylabel('ln|(h(x_{max},t)-h0)/E*h0)|');

figure(3)
plot(time,hmax);
title('x_{max} vs t');
xlabel('t');
ylabel('x_{max}');

