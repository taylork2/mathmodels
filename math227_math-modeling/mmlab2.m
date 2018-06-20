N=100;
array_x = 1:N;
r=3.5; %variable 
x_0=0.5;
array_x(1)=x_0;

for i=1:N
    array_x(i+1) = r*array_x(i)*(1-array_x(i));
%     disp(array_x(i));
end


X=0:0.001:1;
plot(X,X)
hold on
plot(X, r.*X.*(1-X))

for j = 1:N
    plot([array_x(j) array_x(j)], [array_x(j) array_x(j+1)])
    plot([array_x(j) array_x(j+1)], [array_x(j+1) array_x(j+1)])
end

xlabel('x_n')
ylabel('x_{n+1}')
title('r=3.5')
% axis([0 1 0 1])
hold off