clear all;
close all;

h = [.001; .01; .1]; %step size 
e_for = zeros(length(h),1); %error for forward difference scheme 
e_cen = zeros(length(h),1); %error for centered difference scheme 

for i=1:length(h)
    x = 0:h(i):1; 
    y = 2*(cos(pi*x)).^2 - (cos(pi*x)).^4; %y-values 
    dy_exact = -4*pi*sin(pi*x).*cos(pi*x).*(1-(cos(pi*x)).^2); %exact solution for first derivative 

    %Create arrays for numerical sol'ns
    dy_for = x; %forward difference scheme 
    dy_cen = x; %centered difference scheme 
    
    for j=1:(length(x))
        %due to periodic conditions, values "wrap" around 
        if j==length(x)
            dy_for(j) = (y(2) - y(j))/h(i);
            dy_cen(j) = (y(2) - y(j-1))/(2*h(i));
        elseif j==1
            dy_cen(j) = (y(j+1) - y(length(x)-1))/(2*h(i));        
        else 
            dy_for(j) = (y(j+1) - y(j))/h(i);
            dy_cen(j) = (y(j+1) - y(j-1))/(2*h(i));
        end
        
        %Calculating the sigma part of the error 
        e_for(i) = e_for(i)+abs(dy_for(j) - dy_exact(j))^2;
        e_cen(i) = e_cen(i)+abs(dy_cen(j) - dy_exact(j))^2;
    end    
    
    %error calculation 
    N = 1/h(i);
    e_for(i) = sqrt(e_for(i)/N);
    e_cen(i) = sqrt(e_cen(i)/N);

    %Plotting the forward scheme 
    figure
    hold on
    plot(x, dy_exact, 'linewidth', 2);
    plot(x, dy_for, '--', 'linewidth', 2);
    title(['dy/dx with Forward difference scheme where h=', num2str(h(i))]);
    legend('exact', 'numerical', 'location', 'southeast');
    xlabel('x');
    ylabel('dy/dx'); 
    
    %Plotting the centered scheme 
    figure
    hold on
    plot(x, dy_exact, 'linewidth', 2);
    plot(x, dy_cen, '--', 'linewidth', 2);
    title(['dy/dx with Centered difference scheme where h=', num2str(h(i))]);
    legend('exact', 'numerical','location','southeast');
    xlabel('x');
    ylabel('dy/dx'); 
end

%Plotting the error 
figure 
plot(log(h), log(e_for));
title('Forward scheme Error vs h');
xlabel('log(h)');
ylabel('log(error)'); 
coefs = polyfit(log(h), log(e_for), 1); %getting the slope of error 
legend(['slope=',num2str(coefs(1))], 'location','southeast');
figure
plot(log(h), log(e_cen));
title('Centered scheme Error vs h');
xlabel('log(h)');
ylabel('log(error)'); 
coefs = polyfit(log(h), log(e_cen), 1); %getting the slope of error vs h 
legend(['slope=',num2str(coefs(1))], 'location','southeast');