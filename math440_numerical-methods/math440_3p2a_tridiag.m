%2ND CENTRAL DIFFERENCE APPROX TRIDIAG 
clear all;
close all; 

h = [0.025 0.01 0.005]; %step size 

lambda = -1; 
L = 1; %length 

%nonhomogenous part of ode 
f = @(x) x; 

%Exact solution 
y = @(x) csc(1)*sin(x)-x;

%iterate for each step size h 
for i=1:length(h)
    
    x = 0:h(i):L; 
    N = length(x)-2; %number of interior points in x
    u = x; 
    
    %Initial conditions 
    u(1) = 0;
    u(N+2) = 0;

    %Create tridiag matrix 
    A = (1/h(i)^2).*(diag(-1*ones(N-1,1), -1) + diag((2+lambda*h(i)^2)*ones(N,1), 0) + diag(-1*ones(N-1,1), 1));

    %Create f matrix 
    F = f(x(2:N+1));
    F(1) = F(1) + u(1)/h(i)^2; 
    F(N) = F(N) + u(N+2)/h(i)^2;
    
    %linearlly solve the matrices system 
    u(2:N+1) = linsolve(A, F');

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
    err = abs(y(x) - u);
    plot(x, err);
    xlabel('x');
    ylabel('abs error');
    title('Error');
end
