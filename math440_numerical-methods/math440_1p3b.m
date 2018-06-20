clear all; 

h = [.1; 0.05; 0.025; .0125]; %step sizes 

%initialize column variables with same lenghth as h 
un = h;
global_err = h; 
err = h; 

%real solution 
yn = @(t) 1/(t+1);
n=1; %t=1 

for i=1:length(h)
    t = 0:h(i):1; 
    u = t; %initialize u 
    u(1) = 1; %initial condition 
    
    %Euler's method 
    for j = 1:(length(t)-1)
        u(j+1) = u(j) - h(i)*u(j)^2;
    end
    
    un(i) = u(length(t));
    global_err(i) = un(i) - yn(n);
    err(i) = abs(un(i) - yn(n))/h(i); 
    
end

%Create and display the table 
T = table(h, un, global_err, err, 'VariableNames', {'h', 'u_n', 'global_error', 'error'});
disp(T);