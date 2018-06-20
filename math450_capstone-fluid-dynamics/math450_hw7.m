close all;
clearvars;

%initial variables 
N=100; %x steps 
m=50; %number of time steps 
T=1500; %max Time iterations 
M=5000; %max Newton's method iterations  
h0=1; %h initial
E=0.1; %episilon - amplitude of wave 

p=1; %density
g=1; %gravity 
u=100; 
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
h = h_plot(2:N);

hn_plot=h_plot; %create hn 
hn=h;

oldh=h;


plot(x, hn_plot);
hold on

m=10; %# of iterations 

error = zeros(N-1,1);
conv = ones(N-1,1) * 10^-10; %value to determine convergence  
maxval = 10^10; %value to determine divergence 

%Jacobian stuff h(1) = h_i+1 h(2)=h_i h(3) = h_i-1
s=alpha*beta;
df1 = @(h1, h2, h3) s*( 3*(h2+h1)^2 * (h2-h1) - (h2+h1)^3 ); 
df2 = @(h1, h2, h3) 1 - s*( 3*( h3+h2 )^2*( h3-h2 ) - (h3+h2)^3 - ...
    3*(h2+h1)^2*(h2-h1) - (h2+h1)^3 );
df3 = @(h1, h2, h3) -s*( 3*(h3+h2)^2 * (h3-h2) + (h3+h2)^3 ); 
J=zeros(N-1,N-1);

f = @(x, h1, h2, h3) h2 - x - s * ( (h3+h2)^3 * (h3-h2) - (h2+h1)^3*(h2-h1));
b = h; %first guess

    
for t=1:T
    
    for i=1:M
        
        b(1)=f(oldh(1),h0,h(1),h(2));
        b(N-1) = f(oldh(N-1),h(N-2),h(N-1),h0);
        for j=2:N-2
            b(j) = f(oldh(j),h(j-1),h(j),h(j+1));
        end
      
        %Create Jacobian 
        J(1,1) = df2(h0, h(1), h(2)); 
        J(1,2) = df3(h0, h(1), h(2));
        J(N-1,N-1) = df2(h(N-2), h(N-1), h0);
        J(N-1,N-2) = df1(h(N-2), h(N-1), h0);
        for j=2:N-2
            J(j,j-1) = df1(h(j-1),h(j),h(j+1));
            J(j,j) = df2(h(j-1),h(j),h(j+1));
            J(j,j+1) = df3(h(j-1),h(j),h(j+1));
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
    if mod(t,10)==0
        hn_plot(2:N) = h;
        plot(x,hn_plot);
    end
    
end

%graphing the wave 
figure(1)
title('Implicit solution for Height of water');
xlabel('x');
ylabel('height');
