%LAPLACE EQUATION BY FINITE DIFFERENCE

clear all;
close all;

%Analytic Solution 
v = 0;
for m=1:2:801 %(infinity) 
    v = v + 4*sin(m*pi/4)*sinh(m*pi/4)/(m*pi*sinh(m*pi));
end
disp(v);

%initialize h values 
h = [1/4 1/8 1/16 1/32 1/64];
error = zeros(length(h),1);

%iterate over different h values 
for i=1:length(h)
    %Calculate N 
    N = 1/h(i) - 1;
    
    %Create the matrix to plot with initial Dirichlet B.C.s 
    uplot = zeros(N+2, N+2);
    uplot(N+2, :) = 1; %initial condition
    
    %Create b matrix for when v(x,1)=1
    b = zeros(N,N);
    b(N,:) = -1/h(i)^2; 
    b = reshape(b, [], 1);
    
    %initialize matrices to build A matrix
    T = diag(4*ones(N,1)) + diag(-1*ones(N-1,1), -1) + diag(-1*ones(N-1,1), 1);
    
    %Create A matrix 
    A = T;  
    for j=1:N-1
        A = blkdiag(A, T); %This will make T along the whole diag 
    end
    %Add -1's on the Nth and -Nth diagonals
    A = (-1/h(i)^2)*(A + diag(-1*ones(N*N-N,1),-N) + diag(-1*ones(N*N-N,1),N));
    
    %calculate u
    u = A\b;
    
    %Change uplot with calculated u values 
    uplot(2:N+1,2:N+1) = reshape(u, [N,N]);
    
%     %plot 
%     figure 
%     [X,Y] = meshgrid(0:h(i):1); %Create the x and y boundaries 
%     surf(X,Y,uplot);
%     title(['Laplace Equation Surface Plot where h=', num2str(h(i))]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('v');
%     
%     figure 
%     contour(X,Y,uplot);
%     title(['Laplace Equation Contour Plot where h=', num2str(h(i))]);
    
    %calculate the error 
    error(i) = abs(v - uplot((N+1)/4+1,(N+1)/4+1));
end
%plot the error 
figure 
hold on
plot(h, error, '--','linewidth', 3);
title('Error at (1/4, 1/4)');
xlabel('h');
ylabel('error');

coefs = polyfit(h, error', 1);
plot(h, polyval(coefs,h));

coefs = polyfit(h, error', 2);
plot(h, polyval(coefs,h));

legend({'error', 'polyfit 1', 'polyfit 2', 'polyfit 3'}, 'location', 'northwest');

