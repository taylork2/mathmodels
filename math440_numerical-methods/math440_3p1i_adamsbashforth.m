%ADAMS BASHFORTH & ADAMS MOULTON    METHOD 
clear all;
close all; 

%If moulton var. set to true, then will use moulton corrector
%will not otherwise
moulton = true; 

h = [.0125 .025 .05 .1 .2]; %step size 

%First order system equations 
f1 = @(y1,y2) y2;                   %for y
f2 = @(y1,y2) -(y1^2-1)*y2 - y1;    %for y' 

err_max = h; 

%Iterate over different step sizes 
for i = 1:length(h)
    t = 0:h(i):30; %Numerical solution interval 
    
    %initialize y and IC's
    y1 = t; %=y
    y2 = t; %=y' 
    y1(1) = 0.1;
    y2(1) = 0.1; 
    
    %Eulers method for first step 
    y1(2) = y1(1) + h(i)*f1(y1(1),y2(1));
    y2(2) = y2(1) + h(i)*f2(y1(1),y2(1));
    
    %Adams Bashforth for following steps 
    for j=3:length(t)
        %Predictor 2-step AB
        y1(j) = y1(j-1) + (3/2)*h(i)*f1(y1(j-1),y2(j-1)) - (1/2)*h(i)*f1(y1(j-2),y2(j-2));
        y2(j) = y2(j-1) + (3/2)*h(i)*f2(y1(j-1),y2(j-1)) - (1/2)*h(i)*f2(y1(j-2),y2(j-2));
        
        %The corrector if running with 1-step AM 
        if moulton
            y1(j) = y1(j-1) + h(i)/2*(f1(y1(j),y2(j))+f1(y1(j-1),y2(j-1)));
            y2(j) = y2(j-1) + h(i)/2*(f2(y1(j),y2(j))+f2(y1(j-1),y2(j-1)));
        end
    end
    
    %"Real" solution using Matlab ode solver with Runge Kutta  
    [T,Y] = ode45(@vdp1,t,[0.1; 0.1]);
    
    %Plot in physical space 
    figure
    hold on
    plot(T,Y(:,1),'--','linewidth',2);
    plot(t, y1);
    title(['y vs t for h=',num2str(h(i))]);
    xlabel('t');
    ylabel('y');
    legend({'Matlab solution','Numerical solution'}, 'Location','northwest');
    
    %Plot in phase plane 
    figure 
    hold on
    plot(Y(:,1),Y(:,2),'--','linewidth',2);
    plot(y1,y2); 
    title(['y'' vs y for h=',num2str(h(i))]);
    xlabel('y');
    ylabel('y''');    
    legend({'Matlab solution','Numerical solution'}, 'Location','northwest');
    
    figure(3)
    hold on
    plot(t, abs(Y(:,1)'-y1));
    title('Adams Moulton Error vs Time');
    ylabel('Absolute Error');
    xlabel('time');
    axis([0 30 0 5]);
    
%     N = find(t==15);
%     err = abs(Y(:,1)'-y1);
%     err_max(i) = max(err);
%     disp(max(err));
%     disp(sum(err));
end
% 
% figure 
% loglog(h, err_max, 'linewidth', 3);
% % loglog(h, err_max, 'linewidth', 3);
% 
% coefs = polyfit(log(h), log(err_max), 1);
% % hold on
% % plot(h, polyval(coefs,h));
% 
% % coefs = polyfit(h, err_max, 2);
% % plot(h, polyval(coefs,h));
% % 
% % coefs = polyfit(h, err_max, 3);
% % plot(h, polyval(coefs,h));

function dydt = vdp1(T,Y)
    dydt = [Y(2); (1-Y(1)^2)*Y(2)-Y(1)];
end

