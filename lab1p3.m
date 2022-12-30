x_val = 0:0.01:20;
y_val = sin(x_val).*x_val;
plot(x_val,y_val)
y2_val = x_val/2;
hold on
plot(x_val,y2_val,'r')
xlabel('time')
ylabel('position')
title('plot of time vs position')