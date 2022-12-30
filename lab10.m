close all

N=10001;
theta = 1:N;
theta(1) = 1;
theta(2) = 1;
g = 9.80665;
L=10;
dt = 0.001;

for i=3:N
    theta(i) = 2*theta(i-1)-theta(i-2)-dt^2*g*sin(theta(i-1))/L; 
    if i==4001
        force = -0.0022;
        %force = -0.0024;
        theta(i) = theta(i)+force;
    end
end
x=1:N;
y=1:N;

for i=1:N
    x(i)=L*sin(theta(i));
    y(i)=-L*cos(theta(i)); 
    
    %{
    plot([0 x(i)],[0 y(i)], '-p')
    axis([-10 10 -10 10], 'square')
    title(['t=', num2str(i*dt), 'second'])  
    legend((x(i)))
    pause(0.0001)
    %}
end
%hold off
disp(num2str(x(N)))
disp(num2str(y(N)))




    