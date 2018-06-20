%SHOOTING METHOD 
clear all;
close all;

MAXITER = 30;
iter = 1;
err = 1:MAXITER;

h = [0.1 0.05 0.025];

eps = 10^-5; %tolerance 
lambda = 1:MAXITER;

%initial guesses
lambda(1) = -pi; 
lambda(2) = -pi-0.01;
beta = 1:MAXITER; 
m=1;
yL = 0;
L =1;

%First order system equations 
f1 = @(y1,y2) y2;            %for y
f2 = @(y1,y2,lam) lam*y1;    %for y' 

%secant method 
g = @(a0,a1,b0,b1) a1 + (yL - b1)*(a1-a0)/(b1-b0); 

for i=1:length(h)
    t = 0:h(i):L; %Numerical solution interval 
    %initialize y and IC's
    N = length(t);
    y1 = t; %=y
    y2 = t; %=y' 
    y1(1) = 0;
    y2(1) = 1; 
    
    figure(i)
    plot(t, y2(1)*sin(m*pi*t)/(m*pi), '--','linewidth',3);

    while iter < MAXITER 
        if iter>=3
            lambda(iter) = g(lambda(iter-2), lambda(iter-1), beta(iter-2), beta(iter-1));
        end 
        %Eulers method for first step 
        y1(2) = y1(1) + h(i)*f1(y1(1),y2(1));
        y2(2) = y2(1) + h(i)*f2(y1(1),y2(1),lambda(iter));

        %Adams Bashforth for following steps 
        for j=3:N
            %2-step AB
            y1(j) = y1(j-1) + (3/2)*h(i)*f1(y1(j-1),y2(j-1)) - (1/2)*h(i)*f1(y1(j-2),y2(j-2));
            y2(j) = y2(j-1) + (3/2)*h(i)*f2(y1(j-1),y2(j-1),lambda(iter)) - (1/2)*h(i)*f2(y1(j-2),y2(j-2),lambda(iter));
        end

        beta(iter) = y1(N);
        iter = iter+1;

        %Plot in physical space 
        figure(i)
        hold on
        plot(t, y1);
        title(['y vs t for h=',num2str(h(i))]);
        xlabel('t');
        ylabel('y');

        err(iter) = abs(beta(iter) - yL);
        if err(iter)<eps
           break;
        end
    end
    
    disp(lambda(iter-1));

    iter = 1;
    figure(10)
    hold on
    plot(t, abs(y2(1)*sin(m*pi*t)/(m*pi) - y1));
    title('Global Error');
    xlabel('t');
    ylabel('Error');
    
end 