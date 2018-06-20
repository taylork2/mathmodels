close all;
clearvars;

%INITIAL VARIABLES 
N=128; %x steps 
steps=10000; %number of time steps 1000
MAXITER=500; %max Time iterations 1900
Max_Newton=50; %max Newton's method iterations  
h0=.1; %h initial .01
E=0.001; %episilon - amplitude of wave 
L = 1; %length of x 

dx=L/(N-1); dt=.1; %0.01
x = 0:dx:L;
k = 2*pi; %22
h = (h0 + E*sin(k*x))'; %initial values of h
hn=h;
oldh=h;

%to plot the slope 
[MAX,I] = max(h); %location of max of fxn 
y=zeros(1, MAXITER);
time=0:dt:(dt*(MAXITER-1));

p=1; %density
g=0; %gravity is 0 for this one
u=1;  
o=1; %surface tension GRAVITY 0.1 ST 0.001

%van der waals terms 
theta = pi/2; 
m = 2;
n = 3;
M = (n-m)/((m-1)*(n-1)); 
hstar = h0/10; 

%TESTING CONVERGENCE 
error = zeros(N,1);
conv = ones(N,1) * 10^-7; %value to determine convergence  
maxval = 10^3; %value to determine divergence 
flag = 0; %if set to 1, breaks if divergence 

%Coefficients  
alpha  = dt/(8*dx^2); beta=p*g/(3*u); 
alpha2 = dt/(8*dx^4); beta2=o/(3*u); 
alpha3 = dt/(8*dx^2); beta3=1/(3*u);
s=alpha*beta;
v=alpha2*beta2;
v2=alpha3*beta3; 

%van der waals 
gamma = o*(1-cos(theta))/(M*hstar);
dphi =  @(x1,x2) gamma*(-n*hstar^n/((x2+x1)/2)^(n+1) + m*hstar^m/((x2+x1)/2)^(m+1));

%derivatives wrt to x1,x2,x3,x4 and x5 
df1 = @(x1,x2,x3,x4,x5) v*(x3 + x2)^3;
df2 = @(x1,x2,x3,x4,x5) v2*(dphi(x3,x2)*(x3 + x2)^3 - 3*dphi(x3,x2)*(x3 + x2)^2*(x3 - x2)) - v*(3*(x3 + x2)^3 - 3*(x3 + x2)^2*(3*x3 - 3*x2 + x1 - x4) + (x3 + x4)^3);
df3 = @(x1,x2,x3,x4,x5) v*(3*(x3 + x2)^2*(3*x3 - 3*x2 + x1 - x4) + 3*(x3 + x4)^2*(3*x3 - x2 - 3*x4 + x5) + 3*(x3 + x2)^3 + 3*(x3 + x4)^3) - v2*(dphi(x3,x2)*(x3 + x2)^3 + dphi(x4,x3)*(x3 + x4)^3 + 3*dphi(x3,x2)*(x3 + x2)^2*(x3 - x2) + 3*dphi(x4,x3)*(x3 + x4)^2*(x3 - x4)) + 1;
df4 = @(x1,x2,x3,x4,x5) v2*(dphi(x4,x3)*(x3 + x4)^3 - 3*dphi(x4,x3)*(x3 + x4)^2*(x3 - x4)) - v*((x3 + x2)^3 - 3*(x3 + x4)^2*(3*x3 - x2 - 3*x4 + x5) + 3*(x3 + x4)^3);
df5 = @(x1,x2,x3,x4,x5) v*(x3 + x4)^3;

%Gravity scheme (original)
g1 = @(x1,x2,x3,x4,x5) ((x4+x3)^3) *(x4-x3) - ((x3+x2)^3)*(x3-x2);
%Surface Tension scheme
g2 = @(x1,x2,x3,x4,x5) ((x3+x4)^3*(x5-3*x4+3*x3-x2)-(x3+x2)^3*(x4-3*x3+3*x2-x1));
%Van der waals scheme 
g3 = @(x1,x2,x3,x4,x5) ((x3+x4)^3*dphi(x4,x3)*(x4-x3)-(x3+x2)^3*dphi(x3,x2)*(x3-x2));

%Total scheme 
f = @(oldx,x1,x2,x3,x4,x5) -oldx + x3 + s*g1(x1,x2,x3,x4,x5) + v*g2(x1,x2,x3,x4,x5)+ v2*g3(x1,x2,x3,x4,x5);

%Linear stability - Display omega
w = -h0^3/(3*u)*(p*g*k^2 + o*k^4 - k^2*dphi(h0,h0))

%Jacobian 
J=zeros(N,N);
b = h; %first guess

%first plot 
plot(x, h);
hold on

for t=1:MAXITER
    
    for i=1:Max_Newton
        %h(N) and h(N-1) wrap around bc ring topology 
        b(1)=f(oldh(1),h(N-2),h(N-1),h(1),h(2),h(3));
        b(2)=f(oldh(2),h(N-1),h(1),h(2),h(3),h(4));
        %h(1) 
        b(N)   = f(oldh(N),  h(N-2),h(N-1),h(N)  ,h(2),h(3));
        b(N-1) = f(oldh(N-1),h(N-3),h(N-2),h(N-1),h(N),h(2));

        for j=3:N-2
            b(j) = f(oldh(j),h(j-2),h(j-1),h(j),h(j+1),h(j+2));
        end
      
        %Create Jacobian 
        J(1,1) = df3(h(N-2),h(N-1), h(1), h(2), h(3)); 
        J(1,2) = df4(h(N-2),h(N-1), h(1), h(2), h(3));
        J(1,3) = df5(h(N-2),h(N-1), h(1), h(2), h(3));
        J(1,N-2) = df1(h(N-2),h(N-1), h(1), h(2), h(3));
        J(1,N-1) = df2(h(N-2),h(N-1), h(1), h(2), h(3));
        
        J(2,1) = df2(h(N-1), h(1), h(2), h(3), h(4));
        J(2,2) = df3(h(N-1), h(1), h(2), h(3), h(4));
        J(2,3) = df4(h(N-1), h(1), h(2), h(3), h(4));
        J(2,4) = df5(h(N-1), h(1), h(2), h(3), h(4));
        J(2,N-1) = df1(h(N-1), h(1), h(2), h(3), h(4));

        J(N-1,N-3) = df1(h(N-3),h(N-2), h(N-1), h(N), h(2));
        J(N-1,N-2) = df2(h(N-3),h(N-2), h(N-1), h(N), h(2));
        J(N-1,N-1) = df3(h(N-3),h(N-2), h(N-1), h(N), h(2));
        J(N-1,N)   = df4(h(N-3),h(N-2), h(N-1), h(N), h(2));
        J(N-1,2)   = df5(h(N-3),h(N-2), h(N-1), h(N), h(2));

        J(N,N-2) = df1(h(N-2), h(N-1), h(N), h(2), h(3));
        J(N,N-1) = df2(h(N-2), h(N-1), h(N), h(2), h(3));
        J(N,N)   = df3(h(N-2), h(N-1), h(N), h(2), h(3));
        J(N,2)   = df4(h(N-2), h(N-1), h(N), h(2), h(3));
        J(N,3)   = df5(h(N-2), h(N-1), h(N), h(2), h(3));
        
        for j=3:N-2
            J(j,j-2) = df1(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j-1) = df2(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j  ) = df3(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j+1) = df4(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
            J(j,j+2) = df5(h(j-2),h(j-1),h(j),h(j+1),h(j+2));
        end
        
        hn = h - J\b; 
        error = abs(hn-h);
        if error < conv
            h=hn;
            break;
        elseif abs(hn) > maxval
            flag = 1;
            disp("diverge");
            break;
        end     
        h=hn;
    end

    oldh = h;
    
    if mod(t,100) == 0
        plot(x,h);    
    end
    
    %If divergent, then break stop iterating 
    if flag==1
        disp(t);
        break;
    end

    %to get the slope 
    y(t)=log(abs((h(I)-h0)/(E*h0)));
end

%graphing the wave 
figure(1)
title('Implicit Height of water w/ van der waals');
xlabel('x');
ylabel('height');


figure(2)
plot(time, y);
title('Slope at x_{max} vs t');
xlabel('t');
ylabel('y');

coef = polyfit(time, y, 1);
disp(coef(1));