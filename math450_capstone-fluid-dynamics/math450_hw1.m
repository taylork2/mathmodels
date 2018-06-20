close all;

N = 8;
m = 100; 
dt = 0.001;
dx = 1/N; 

x=linspace(0,1,N+2); %Create the x array 

T = ones(1,N+2,1); %Create array of initial values which are all 1
%T = zeros(1,N+2,1); %Create array of inital values which are all 0

T(N+2) = 0; %Initially last value is 0
%T(1) = 1; %Initially first value is 1

T_n = T; %Initialize T_n

for j=1:m %number of iterations 
    
    for i=2:N+1 
        T_n(i) = T(i) + dt/dx^2 * (T(i-1) - 2*T(i) + T(i+1));   
    end
    
    if mod(j,25)==0
        plot(x, T_n(1:length(T_n)));
        hold on;    
    end
    T=T_n; %Set new values of T_n to T for next iteration 

end

xlabel('x');
ylabel('f(x)');
title('Numerical Solution');