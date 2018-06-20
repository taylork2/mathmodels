%%
clear all; 
close all;

h0 = 10*10^-9; %[m]
hstar = [h0/10, h0/100, h0/1000];

N=10^5; 

m = [2 3 3];
n = [3 4 9];

h=linspace(0, 10^-8, N);

for i=1:length(hstar) 
    for j=1:length(m)
        phi = (hstar(i)./h).^n(j) - (hstar(i)./h).^m(j); 
        dphi = (-n(j)*hstar(i)^n(j))./(h.^(n(j)+1)) + (m(j)*hstar(i)^m(j))./(h.^(m(j)+1));
        figure
        subplot(2,1,1);
        plot(h, phi);
        title(['Phi vs h where (m,n)=',num2str(m(j)),',',num2str(n(j)),' and hstar=', num2str(hstar(i))]);
        %constrain boundaries to better display "bump"
        if i==3
            xlim([0 10^-10]);
        elseif i==2
            xlim([0 10^-9]);    
        end 
        ylim([-inf 1]); %set the ymin to the automatically calculated min 
        xlabel('h');
        ylabel('Phi');
        subplot(2,1,2);
        plot(h, dphi);
        title(['dPhi/dh vs h where (m,n)=',num2str(m(j)),',',num2str(n(j)),' and hstar=', num2str(hstar(i))]);
        %constrain boundaries to better display "bump"
        if i==3
            xlim([0 10^-10]);
        elseif i==2
            xlim([0 10^-9]);
        end 
        ylim([-1 inf]); %set the ymax to the automatically calculated max 
        xlabel('h');
        ylabel('dPhi');
    end
end

%%
clear all;
close all;
p = 1; %
g = 0; %gravity 
mu = 1;  
sigma = 1;
theta = pi/180.*[15 30 45]; %changing theta from 15-45deg 
m = [2 3 3];
n = [3 4 9];
M = @(n, m) (n-m)/((m-1)*(n-1));

h0 = 0.01;
hstar = [h0/10 h0/100 h0/1000]; %change this 
eps = .001; 
N = 1000;

k = linspace(0, 80, N);

% alpha = p*g*k.^2;
% beta = sigma * k.^4; 
% gamma = k.^2*sigma * (1 - cos(theta))/(M*hstar); 
% dPhi = -n * hstar^n/h0^(n+1) + m*hstar^m/h0^(m+1);
% 
% omega = -h0^3/(3*mu)*(alpha + beta - gamma*dPhi); 

omega = @(theta, m, n, hstar) -h0^3/(3.*mu)*((p*g*k.^2) + sigma*k.^4 - k.^2*sigma*(1 - cos(theta))/(M(n, m)*hstar)*(-n * hstar^n/h0^(n+1) + m*hstar^m/h0^(m+1)));

figure(1)
hold on
for i=1:length(theta) 
    o = omega(theta(i), m(1), n(1), hstar(1));
    plot(k, o);
end
title('w vs k with changing theta');
legend('15deg', '30deg', '45deg');
xlabel('k');
ylabel('w');
ylim([-0.5 0.5]);

figure(2)
hold on
for i=1:length(m)
    o = omega(theta(1), m(i), n(i), hstar(1));
    plot(k, o);
end
title('w vs k with changing (n,m)');
legend('(m,n)=(2,3)', '(m,n)=(3,4)', '(m,n)=(3,9)');
xlabel('k');
ylabel('w');
ylim([-0.01 0.01]);

figure(3)
hold on 
for i=1:length(hstar) 
    o = omega(theta(1), m(1), n(1), hstar(i));
    plot(k, o);
end
title('w vs k with changing hstar');
legend('h0/10', 'h0/100', 'h0/1000');
xlabel('k');
ylabel('w');
ylim([-0.01 0.01]);


