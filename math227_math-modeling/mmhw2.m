x = [1.0 2.3 3.7 4.2 6.1 7.0];
y = [3.6 3.0 3.2 5.1 5.3 6.8];

xb = [29.1 48.2 72.7 92.0 118 140 165];
yb = [0.0493 0.0821 0.123 0.154 0.197 0.234 0.274];

%least squares 

x1=linspace(0,10);
plot(x,y, 'o');
hold on;
p=polyfit(x,y,1);
f1= polyval(p,x1);
plot(x1,f1,'-')


n = 1;
m = length(x);

z = ((x-min(x))-(max(x)-x))/(max(x)-min(x));

A = ones(1,m);
if n > 1
   A(:,2) = z;
end
if n > 2
  for k = 3:n+1
     A(:,k) = 2*z.*A(:,k-1) - A(:,k-2);  %% recurrence relation
  end
end

b=A./y


p2=polyfit(x, b, 1);
f2= polyval(p2,x1);
plot(x1,f2,'-')
