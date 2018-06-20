%FTCS
clear all;
close all;

dx = 0.05;
dt = [0.0013, 0.0012];

t = [0, 1, 25, 50]; %times to plot at 

x = 0:dx:1;
N = length(x)-2;
v = zeros(length(x),max(t)); %initialize v 
v(:,1) = 1-2*abs(x-1/2); %initial condition 

legendInfo = cell(length(t),1);

%Create the diag(1, -2, 1) matrix 
A = (1/dx^2).*(diag(1*ones(N-1,1), -1) + diag(-2*ones(N,1), 0) + diag(1*ones(N-1,1), 1));

%iterate over different dt 
for i=1:length(dt)  
    figure
    n=1;%Index for labeling the legend 
    
    %time iterations 
    for k=1:max(t)+1
        
        %FTCS method 
        v(2:N+1,k+1) = v(2:N+1,k) + dt(i).*A*v(2:N+1,k);
        
        %Plot only if in t 
        if ismember(k-1,t)
            hold on
            plot(x,v(:,k), '-x'); %Numerical solution 
            legendInfo{n} = ['t=',num2str(t(n)),'dt'];
            n=n+1;
        end
        
    end
    %Labeling the graph 
    legend(legendInfo);
    title(['Heat eqn FTCS dt=', num2str(dt(i))]);
    xlabel('x');
    ylabel('v');
end
