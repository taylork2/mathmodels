%NEUMANN TRIDIAG USING GHOST PTS
clear all;
close all; 

h = [0.1 0.01 0.005];

lambda = -1;
L = 1;

%nonhomogenous part of ode 
f = @(x) x;

%Exact solution 
y = @(x) sin(x)-tan(.5)*cos(x)-x;

for i=1:length(h)
    x = 0:h(i):L; 
    N = length(x); %number of interior points in x

    %Create tridiag matrix 
    A = (1/h(i)^2).*(diag(-1*ones(N-1,1), -1) + diag((2+lambda*h(i)^2)*ones(N,1), 0) + diag(-1*ones(N-1,1), 1));
    A(1,2) = -2/h(i)^2; 
    A(N,N-1) = -2/h(i)^2; 

    %Create f matrix 
    F = f(x);

    %linearlly solve the matrices system 
    u = linsolve(A, F');

    %plot the numerical solution 
    figure 
    hold on
    plot(x, y(x), 'linewidth', 1);
    plot(x, u,'--','linewidth', 3);
    xlabel('x');
    ylabel('y');
    legend({'Exact solution','Numerical solution'}, 'location','northwest');
    title(['y vs x with h=',num2str(h(i))]);
    
    figure(2)
    hold on 
    err = abs(y(x) - u');
    plot(x, err);
    xlabel('x');
    ylabel('abs error');
    title('Error');
end
