close all;
clearvars;

N = 5;
k1=1;
k2=2;

%construct the diag matrix 
a1 = ones(1, N^2-1, 1);
for i=1:N^2-1
    if mod(i, N) == 0
        a1(1, i) = 0;
    end
end
a2 = ones(1, N^2-N, 1);
A = diag(a2, N) + diag(-4*ones(1, N^2, 1)) + diag(a2, -N) + diag(a1, 1) + diag(a1, -1);

%initial conditions 

change=floor(N^2/2);
T = zeros(N^2, 1);
for j=1:N:N^2
    T(j)=-k1;
end

%calulate the T values 
T_n = A\T;

%Create matrix including boundary conditions for graphing 
T_plot = zeros(N+2);
T_plot(1:N+2, 1) = 1;
T_n = reshape(T_n,[N,N]);
T_plot(2:N+1, 2:N+1) = rot90(T_n,1); 

%Graph the results 
x = linspace(0,1,N+2); 
y = linspace(0,1,N+2);
[X, Y] = meshgrid(x,y);
surf(X,Y,T_plot);
    
title('2d Laplace(N=100)');

