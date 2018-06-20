close all;
clearvars;

N = 100;
m = 200; 
dt = 0.01;
dx = 1/N; 
k1=1;
k2=2;
shift = N/2; %Find where the shift will occur
ushift = ceil(shift); %upper shift, if N is odd
lshift = floor(shift); %lower shift, if N is odd

x = linspace(0,1,N+1); %Create the x array 

A = zeros(N-1); %Create array of values initially which are 0
T = ones(N-1,1,1);
T(N-1) = 0;
T_n_plot = ones(N+1, 1, 1);
T_n_plot(N+1) = 0;

%Create diag matrix 
for i=1:N-1
    for j=1:N-1   
        if j<lshift
            k=k1;
            K=k1;
        elseif j>ushift
            k=k2;
            K=k2;
        elseif i<shift
            K=k1;
        else
            k=(k1+k2)/2;
            K=k2;
        end
        
        if i==j
            A(i,j) = 1 + 2*k*dt/dx^2;
        elseif i==j+1 || j==i+1
            A(i,j) = -K*dt/dx^2;
        end
    end
end

%Calculate the new values 
for j=1:m
    T(1) = T(1)+dt/dx^2;
    T_n = linsolve(A,T);
    T_n_plot(2:N) = T_n(1:N-1);
    plot(x, T_n_plot);
    hold on
    T=T_n;
end
    
xlabel('x');
ylabel('f(x)');
title('Matrix Numerical Solution, k1=1, k2=2');
