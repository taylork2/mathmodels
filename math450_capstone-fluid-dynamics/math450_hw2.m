close all;
clearvars;

N = 50;
m = 500; 
dt = 0.01;
dx = 1/N; 
k=2;
shift = floor(N/2);

x = linspace(0,1,N+1); %Create the x array 

A = zeros(N-1); %Create array of values initially which are 0
A2 = zeros(N-1);
T = ones(N-1,1,1);
T2 = ones(N-1,1,1);
T2(N-1) = 0;
T(N-1) = 0;
T_n_plot = ones(N+1, 1, 1);
T_n_plot(N+1) = 0;

%Create diag matrix 
for i=1:N-1
    for j=1:N-1
        if i==j
            A(i,j) = 1 + 2*dt/dx^2;
        elseif i==j+1 || j==i+1
            A(i,j) = -dt/dx^2;
        end
    end
end

for i=1:N-1
    for j=1:N-1
        if i==j
            A2(i,j) = 1 + 2*k*dt/dx^2;
        elseif i==j+1 || j==i+1
            A2(i,j) = -k*dt/dx^2;
        end
    end
end

for j=1:m
    T(1) = T(1)+dt/dx^2;
    T2(1) = T2(1)+k*dt/dx^2;
    %T(N-1) = dt/dx^2;
    T_n = linsolve(A,T);
    T_n2 = linsolve(A2,T2);
    %T_n_plot(2:shift) = T_n(1:shift-1);
    %T_n_plot(shift+1:N) = T_n2(shift:N-1);
    plot(x, T_n_plot);
    hold on
    T=T_n;
end
    
xlabel('x');
ylabel('f(x)');
title('Matrix Numerical Solution');
