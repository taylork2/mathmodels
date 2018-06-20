close all;
clearvars;

%initial variables 
N=100; %x steps 
m=10; %number of time steps 
T=2500; %max Time iterations 
M=50; %max Newton's method iterations  
h0=1; %h initial 
E=0.1; %episilon - amplitude of wave 

p=1; %density
g=1; %gravity 
u=1000; 
k=2*pi; 
o=0.00001; %surface tension
beta=p*g/(3*u); 
delta=o/(3*u); 

dx=1/N;
dt=1/m;
alpha = dt/(8*dx^2); 
gamma = dt/(8*dx^4);

x = linspace(0,1,N+1); %initialize x variables
h_plot = (h0 + E*sin(k*x))'; %initial values of h
h = h_plot(2:N);

hn_plot=h_plot; %create hn 
hn=h;

oldh=h;

plot(x, hn_plot);
hold on

m=10; %# of iterations 

error = zeros(N-1,1);
conv = ones(N-1,1) * 10^-7; %value to determine convergence  
maxval = 10^3; %value to determine divergence 

%Jacobian stuff h(1) = h_i+1 h(2)=h_i h(3) = h_i-1
s=alpha*beta;
v=gamma*delta;
df1 = @(h1, h2, h3, h4, h5) -v*(h3+h2)^3; 

df2 = @(h1, h2, h3, h4, h5) s*( 3*(h3+h2)^2 * (h3-h2) - (h3+h2)^3 ) + ...
    v*( 3*(h3+h2)^3 + 3*(h3+h2)^2*(h4-3*h3+3*h2-h1) + (h4+h3)^3); 

df3 = @(h1, h2, h3, h4, h5) 1 - s*( 3*( h4+h3 )^2*( h4-h3 ) - (h4+h3)^3 - ...
    3*(h3+h2)^2*(h3-h2) - (h3+h2)^3 ) - v*( 3*(h4+h3)^3 + 3*(h4+h3)^2*(h5-3*h4+3*h3-h2) - ... 
    ( -3*(h3+h2)^3 + 3*(h3+h2)^2*(h4-3*h3+3*h2-h1) ));

df4 = @(h1, h2, h3, h4, h5) -s*( 3*(h4+h3)^2 * (h4-h3) + (h4+h3)^3 ) - ...
    v * ( -3*(h4+h3)^3 + 3*(h4+h3)^2*(h5-3*h4+3*h3-h2) - (h3+h2)^3); 

df5 = @(h1, h2, h3, h4, h5) -v * (h4+h3)^3; 
J=zeros(N-1,N-1);

f = @(x, h1, h2, h3, h4, h5) h3 - x - s * ( (h4+h3)^3 * (h4-h3) - (h3+h2)^3*(h3-h2)) - ...
    v * ((h4+h3)^3*(h5- 3*h4 + 3*h3 -h2) - (h3+h2)^3*(h4- 3*h3 + 3*h2 -h1));
b = h; %first guess

    
for t=1:T
    
    for i=1:M
        %h(N-1) is ghost point 
        b(1)=f(oldh(1),h(N-1),h0,h(1),h(2),h(3));
        b(2)=f(oldh(2),h0,h(1),h(2),h(3),h(4));
        %h(1) is ghost point 
        b(N-1) = f(oldh(N-1),h(N-3),h(N-2),h(N-1),h0,h(1));
        b(N-2) = f(oldh(N-2),h(N-4),h(N-3),h(N-2),h(N-1),h0);

        for j=3:N-3
            b(j) = f(oldh(j),h(j-2),h(j-1),h(j),h(j+1),h(j+2));
        end
      
        %Create Jacobian 
        J(1,1) = df3(h(N-1),h0, h(1), h(2), h(3)); 
        J(1,2) = df4(h(N-1),h0, h(1), h(2), h(3));
        J(1,3) = df5(h(N-1),h0, h(1), h(2), h(3));

        J(2,1) = df2(h0, h(1), h(2), h(3), h(4));
        J(2,2) = df3(h0, h(1), h(2), h(3), h(4));
        J(2,3) = df4(h0, h(1), h(2), h(3), h(4));
        J(2,4) = df5(h0, h(1), h(2), h(3), h(4));

        J(N-2,N-4) = df1(h(N-4),h(N-3),h(N-2), h(N-1), h0);
        J(N-2,N-3) = df2(h(N-4),h(N-3),h(N-2), h(N-1), h0);
        J(N-2,N-2) = df3(h(N-4),h(N-3),h(N-2), h(N-1), h0);
        J(N-2,N-1) = df4(h(N-4),h(N-3),h(N-2), h(N-1), h0);        

        J(N-1,N-3) = df1(h(N-3),h(N-2), h(N-1), h0, h(1));
        J(N-1,N-2) = df2(h(N-3),h(N-2), h(N-1), h0, h(1));
        J(N-1,N-1) = df3(h(N-3),h(N-2), h(N-1), h0, h(1));
        
        for j=3:N-3
            J(j,j-2) = df1(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j-1) = df2(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j  ) = df3(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j+1) = df4(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j+2) = df5(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
        end
        
        hn = h - J\b; 
        error = abs(hn-h);
        if error < conv
            h=hn;
            break;
        elseif abs(hn) > maxval
            break;
        end     
        h=hn;
    end
    
    oldh = h;
    
    %prevent each step from being graphed (for computational time's sake)
    if mod(t,50)==0
        hn_plot(2:N) = h;
        plot(x,hn_plot);
    end
    
end

%graphing the wave 
figure(1)
title('Implicit solution for Height of water with Surface Tension');
xlabel('x');
ylabel('height');
