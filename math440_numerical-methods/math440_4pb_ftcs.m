%FTCS
clear all;
close all;

dx = 0.05;
dt = [0.0012, 0.0013];

t = [0, 1, 25, 50, 550]; %times to plot at 

x = 0:dx:1;
N = length(x)-2;
v = zeros(length(x),max(t)); %initialize v 
v(:,1) = sin(pi*x); %initial condition 

err = zeros(max(t),1); 

legendInfo = cell(length(t)+1,1);
errLegendInfo = cell(length(dt),1);

%Create the diag(1, -2, 1) matrix 
A = (1/dx^2).*(diag(1*ones(N-1,1), -1) + diag(-2*ones(N,1), 0) + diag(1*ones(N-1,1), 1));

%iterate over different dt 
for i=1:length(dt)  
    figure
    n=1; %This index will be for labeling the graph
    
    %time iterations 
    for k=1:max(t)+1
        
        %FTCS method 
        v(2:N+1,k+1) = v(2:N+1,k) + dt(i).*A*v(2:N+1,k);
        
        %Plot only if in t 
        if ismember(k-1,t)
            hold on
            plot(x, exp(-pi^2*(k-1)*dt(i))*sin(pi*x), 'linewidth', 2); %Exact solution 
            plot(x,v(:,k), '-x'); %Numerical solution 
            
            %Creating the legend 
            legendInfo{2*n-1} = ['Exact t=',num2str(t(n)),'dt'];
            legendInfo{2*n} = ['t=',num2str(t(n)),'dt'];
            n=n+1;
        end
        
    end
    %Labeling the graph 
    legend(legendInfo);
    title(['Heat Eqn FTCS dt=', num2str(dt(i))]);
    xlabel('x');
    ylabel('v');
    
    %Plotting the error 
    figure(3)
    hold on
    for j=1:length(err)
        err(j) = sum(abs(v(:,j)' - exp(-pi^2*(j+1)*dt(i))*sin(pi*x)));   
    end
    
    errLegendInfo{i} = ['dt=',num2str(dt(i))];
    plot(0:length(err)-1, err);
    
end

title('Total Absolute Error');
xlabel('t');
ylabel('Total Abs error');
legend(errLegendInfo);
