clear all;
close all;

m=20; %# of iterations 
x0 = [ 1; 1; 1 ]; %initial conditions 

x = x0;
xn = x;
error = zeros(3,1);
conv = ones(3,1) * 10^-10; %value to determine convergence  
maxval = 10^10; %value to determine divergence 

%define functions 
% f1 = @(x) 1-2*x(1)+x(2)+x(1)*((x(2)-1)/2);
% f2 = @(x) x(1)-2*x(2)+x(3)+x(2)*((x(3)-x(1))/2);
% f3 = @(x) x(2)-2*x(3)+x(3)*(-x(2)/2);
f1 = @(x) 3*x(1) - cos(x(2)*x(3)) - 1/2;
f2 = @(x) x(1)^2 - 81*(x(2)+0.1)^2 + sin(x(3)) + 1.06;
f3 = @(x) exp(-x(1)*x(2)) + 20*x(3) + 1/3 * (10*pi - 3);

%create array of function values 
f = @(x) [ f1(x); f2(x); f3(x) ];

%partial derivatives of above functions 
% df11 = @(x) -2 + (x(2)-1)/2;
% df12 = @(x) 1 + x(1)/2;
% df21 = @(x) 1 - x(2)/2;
% df22 = @(x) -2 + (x(3)-x(1))/2;
% df23 = @(x) 1 + x(2)/2;
% df32 = @(x) 1 - x(3)/2;
% df33 = @(x) -2 - x(2)/2;
df11 = @(x) 3;
df12 = @(x) x(3)* sin(x(3)*x(2));
df13 = @(x) x(2) * sin(x(3)*x(2));
df21 = @(x) 2*x(1);
df22 = @(x) -2*81*(x(2)+0.1);
df23 = @(x) cos(x(3));
df31 = @(x) -x(2)*exp(-x(1)*x(2));
df32 = @(x) -x(1)*exp(-x(1)*x(2));
df33 = @(x) 20;

%create inverse Jacobian matrix 
J = @(x) [ df11(x) df12(x) df13(x); df21(x) df22(x) df23(x); df31(x) df32(x) df33(x) ]';

hold on
plot(0:1:2, xn, '-o');

for i=1:m
    disp(['Iteration ', num2str(i)]);
    xn = x - J(x)\f(x)
    error = abs(xn-x);
    if error < conv
        disp('converged');
        break;
    elseif abs(xn) > maxval
        disp('diverged');
        break;
    end
    hold on
    plot(0:1:2, xn, '-o');
    x=xn;
end

axis([0 2 -0.6 1]);
title('Newtons Method on System of Equations');
xticks([0 1 2]);
xlabel('zeroes');
ylabel('x value');


