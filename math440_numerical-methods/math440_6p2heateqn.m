%2-D HEAT EQN OVER TIME 
clear all;
close all;

%initialize values 
h = 0.1;
dt = 0.0025;
ti=0; tf=2; 

%Calculate N, size of matrix 
N = 1/h - 1;

%Create the matrix to plot with initial Dirichlet B.C.s 
uplot = zeros(N+2, N+2);
uplot(N+2, :) = 1; %initial condition

%Create matrix for u at t=0
u0 = zeros(N,N);
u0 = reshape(u0, [], 1);

%Create b matrix 
b = zeros(N,N);
b(N,:) = dt/h^2; 
b = reshape(b, [], 1);

%initialize matrices to build A matrix
T = diag(4*ones(N,1)) + diag(-1*ones(N-1,1), -1) + diag(-1*ones(N-1,1), 1);

%Create A matrix 
A = T;  
for j=1:N-1
    A = blkdiag(A, T); %This will make T along the whole diag 
end
%Add -1's on the Nth and -Nth diagonals
A = (-dt/h^2)*(A + diag(-1*ones(N*N-N,1),-N) + diag(-1*ones(N*N-N,1),N));

%iterate over time 
for t=ti:dt:tf 
    %calculate u
    u = u0 + A*u0 + b;
    u0=u;
end

%Change uplot with calculated u values 
uplot(2:N+1,2:N+1) = reshape(u, [N,N]);

%plot 
figure 
[X,Y] = meshgrid(0:h:1); %Create the x and y boundaries 
surf(X,Y,uplot);
title(['2D Heat Equation Surface Plot where t=', num2str(t)]);
xlabel('x');
ylabel('y');
zlabel('v');

figure 
contour(X,Y,uplot);
title(['2D Heat Equation Contour Plot where t=', num2str(t)]);

