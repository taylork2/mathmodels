%UPWIND/DOWNWIND SCHEME 
clear all;
close all;

syms X

scheme = 'u'; %choose which scheme to use 
F = '2'; %choose which initial data to use 

%initial conditions 
dx = 0.05;
x = -1:dx:5;
dt = [0.04];
ti = 0;
tf = 2;

%initial data at v(x,0)
f1=piecewise(X<0, 1, X==0, 0, X>0, -1);
f2=piecewise(X<0, -1, 0<=X<=2, 1-2*abs(X-1), X>2, -1);

%set f to correct initial data
if F=='1'
    f = f1; 
elseif F=='2'
    f=f2;
end

err = dt;

%iterate over each dt 
for i=1:length(dt) 
    t = ti:dt(i):tf;  %initialize time vector 
    u = zeros(length(x), length(t)); %initialize space vector 
    
    %set initial data 
    u(:,1) = subs(f, X, x);  
%     figure 
%     hold on
%     %plot real solutions 
%     plot(x, subs(f,X,x-ti), 'linewidth', 2);
%     plot(x, subs(f,X,x-tf), 'linewidth', 2);
    
    %iterate over time 
    for n=1:length(t)-1 
        u(1,n+1) = u(1,n); %initial point 
        u(length(x), n+1) = u(length(x), n); %last point 
        %upwind scheme 
%         if scheme=='u'
            for j=2:length(x)
                u(j,n+1) = (1 - dt(i)/dx)*u(j,n) + dt(i)/dx*u(j-1,n);
            end
            uu = u;
        %downwind scheme 
%         elseif scheme=='d'
            for j=1:length(x)-1
                u(j,n+1) = (1 + dt(i)/dx)*u(j,n) - dt(i)/dx*u(j+1,n);
            end
            ud = u;
        %lax-friedrichs     
%         elseif scheme=='f'
            for j=2:length(x)-1
                u(j,n+1) = (1/2)*(u(j+1,n)+u(j-1,n)) - (dt(i)/(2*dx))*(u(j+1,n)-u(j-1,n));
            end
            uf = u;
        %lax-wendroff 
%         elseif scheme=='w'
            for j=2:length(x)-1
                u(j,n+1) = u(j,n) - (dt(i)/(2*dx))*(u(j+1,n)-u(j-1,n)) + (1/2)*(dt(i)/dx)^2*(u(j+1,n)-2*u(j,n)+u(j-1,n));
            end    
            uw = u; 
%         end      
    end
    
    %plot at last time step 
%     hold on
%     plot(x,u(:,n+1), '-o');
    
%     %label plot 
%     xlabel('x');
%     ylabel('u');
%     legend(['Exact t=', num2str(ti)], ['Exact t=', num2str(tf)], 'Num tf');
%     if scheme=='d'
%         title(['Downwind scheme of f', F, '(x) dt=', num2str(dt(i))]);
%     elseif scheme=='u'
%         title(['Upwind scheme of f', F, '(x) dt=', num2str(dt(i))]);    
%     elseif scheme=='f'
%         title(['Lax-Friedrichs scheme of f', F, '(x) dt=', num2str(dt(i))]);    
%     elseif scheme=='w'
%         title(['Lax-Wendroff scheme of f', F, '(x) dt=', num2str(dt(i))]);    
%     end
    
    figure
    hold on 
    er = abs(subs(f,X,x-tf)' - u(:,n+1));
    err(i) = sum(er);
%     plot(x, er);
    %Sum of Total Absolute Error 
    disp(sum(abs(subs(f,X,x-tf)' - uu(:,n+1))));
    disp(sum(abs(subs(f,X,x-tf)' - uf(:,n+1))));
    disp(sum(abs(subs(f,X,x-tf)' - uw(:,n+1))));

    plot(x, abs(subs(f,X,x-tf)' - uu(:,n+1)));
%     plot(x, abs(subs(f,X,x-tf)' - ud(:,n+1)));
    plot(x, abs(subs(f,X,x-tf)' - uf(:,n+1)));
    plot(x, abs(subs(f,X,x-tf)' - uw(:,n+1)));
end