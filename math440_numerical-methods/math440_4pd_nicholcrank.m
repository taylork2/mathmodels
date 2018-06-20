%CRANK NICOLSON
clear all;
close all;

dx = 0.05;
dt = [0.0012, 0.0013];

x = 0:dx:1;
N = length(x)-2;

% Create the tridiag(1, -2, 1) matrix  
A = (1/dx^2).*(diag(1*ones(N-1,1), -1) + diag(-2*ones(N,1), 0) + diag(1*ones(N-1,1), 1));

%iterate over different dt
for i=1:length(dt)
    t = 0:dt(i):50*dt(i); %initialize time 
    v = zeros(length(x),length(t)); %initialize v 
    v(:,1) = 1-2*abs(x-1/2); %initial condition 

    figure    
    plot(x, v(:,1), 'linewidth',2); %plot initial condition 
    
    %Crank nicolson method 
    for k=1:length(t)
        v(2:N+1,k+1) = (eye(N)-dt(i)/2*A)^-1*(eye(N)+dt(i)/2*A)*v(2:N+1,k);
        
        hold on
        plot(x, v(:,k+1)); %plot the numerical solution at each time step 
    end
    
    title(['Crank Nicolson dt=', num2str(dt(i))]);
    xlabel('x');
    ylabel('v');
end