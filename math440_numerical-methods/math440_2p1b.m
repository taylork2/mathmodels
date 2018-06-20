%BACKWARD EULER'S AND RUNGE KUTTA 

clear all; 
close all;

% h = [0.15 0.13 0.11 0.10 0.09 0.06 0.05 0.01]; %step sizes 
h = 0.01:0.001:0.3;

%real solutions 
y1 = @(t) -0.1*exp(-20*t) + 1.1*exp(-2*t);
y2 = @(t) 0.1*exp(-20*t) + 1.1*exp(-2*t);

%initial value equations 
dy1 = @(x1,x2) -11*x1 + 9*x2;
dy2 = @(x1,x2) 9*x2 - 11*x2;

%For backwards Euler's linsolve   
A = [-11 9; 9 -11];
I = eye(2); %identity matrix 

legendInfo = cell(length(h), 1);
err_eul = length(h);
err_run = length(h);

for i=1:length(h)
    t = 0:h(i):2; %time steps 
    u1 = t; %initialize u1
    u2 = t; %initialize u2 

    
    %This is a temp variable to hold array of next point in u1 & u2
    %for linearly solving equation 
    u = [0; 0]; 
    
    %initial conditions 
    u1(1) = 1; 
    u2(1) = 1.2; 
    
    err = t;
    err2 = t;
    
    %Backward Euler's scheme 
    for j=1:length(t)-1
        u = linsolve( (I-h(i)*A), [u1(j); u2(j)]);
        u1(j+1) = u(1);
        u2(j+1) = u(2);
        
        err(j+1) = abs(y1(t(j+1)) - u1(j+1)); 
    end

%     figure
%     hold on
%     plot(t, y1(t),'linewidth',3);
%     plot(t, u1, 'linewidth',2);
%     title(['Backward Euler: y1 vs t where h=',num2str(h(i))]);
%     legend('Exact','Numerical');
%     xlabel('time (t)');
%     ylabel('y1');
%     
    err_eul(i) = max(err);
%     legendInfo{i} = ['h=',num2str(h(i))];

    %runge kutta method 
    for j=1:length(t)-1
        k1_1 = dy1(u1(j),u2(j));
        k1_2 = dy2(u1(j),u2(j));
        
        k2_1 = dy1(u1(j)+h(i)*k1_1/2, u2(j)+h(i)*k1_2/2);
        k2_2 = dy2(u1(j)+h(i)*k1_1/2, u2(j)+h(i)*k1_2/2);
        
        k3_1 = dy1(u1(j)+h(i)*k2_1/2, u2(j)+h(i)*k2_2/2);
        k3_2 = dy2(u1(j)+h(i)*k2_1/2, u2(j)+h(i)*k2_2/2);        
        
        k4_1 = dy1(u1(j)+h(i)*k2_1/2, u2(j)+h(i)*k2_2/2);
        k4_2 = dy2(u1(j)+h(i)*k2_1/2, u2(j)+h(i)*k2_2/2);        
        
        u1(j+1) = u1(j)+h(i)*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
        u2(j+1) = u2(j)+h(i)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;  
        
        err2(j+1) = abs(y1(t(j+1)) - u1(j+1)); 
        
    end
    
%     figure
%     hold on
%     plot(t, y1(t),'linewidth',3);
%     plot(t, u1, 'linewidth',2);
%     title(['Runge Kutta: y1 vs t where h=',num2str(h(i))]);
%     legend('Exact','Numerical');
%     xlabel('time (t)');
%     ylabel('y1');
    
%     figure(20)
%     hold on
%     plot(t, err2); 
    err_run(i) = max(err2);
    
end

figure(21)
plot(h, err_eul);
% legend(legendInfo);
title('Backward Euler Error vs t');
xlabel('time (t)');
ylabel('Absolute error');

figure(20)
plot(h, err_run);
% legend(legendInfo);
title('Runge Kutta Error vs t');
xlabel('time (t)');
ylabel('Absolute error');



