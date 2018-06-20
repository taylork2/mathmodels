%LEAP FROG METHOD 

clear all; 
close all;
 
h = [0.5 .25 .125]; %time step size 
lambda = [1 -1]; %changing lambda 

y = @(t, l) exp(l*t);

legendInfo = cell(length(h)*length(lambda),1);

for i=1:length(h)
    
    t = 0:h(i):5; %time 
    u = t; %initialize u with length of time  
    
    err = t; %initialize error 

    for l = 1:length(lambda)       
        %initial conditions 
        u(1) = 1;
        u(2) = exp(h(i)*lambda(l));
        
        %initial error 
        err(1) = 0;
        err(2) = abs(u(2) - y(t(2),lambda(l)));
        
        %Leap-frog method 
        for j = 2:length(t)-1
            u(j+1) = u(j-1) + 2*h(i)*lambda(l)*u(j);
            err(j+1) = abs(u(j+1) - y(t(j+1),lambda(l))); 
        end
        
        %Plotting the graphs 
        figure 
        hold on
        plot(t, y(t,lambda(l)), 'linewidth', 3);
        plot(t, u,'linewidth',2); 
        title(['Leap-frog: y vs t where lambda=',num2str(lambda(l)),' h=',num2str(h(i))]);
        xlabel('time (t)');
        ylabel('y');
        legend('Exact', 'Numerical');
        
        figure(4)
        hold on
        plot(t, err);
        legendInfo{(i-1)*length(lambda)+l} = ['lambda=',num2str(lambda(l)),' h=',num2str(h(i))];

    end    
end

figure(4)
legend(legendInfo);
title('Absolute error');
xlabel('time (t)');
ylabel('Error');