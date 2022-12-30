x = linspace(-3, 3);

sinf1 = 2/pi*sin(pi*x);
sinf2 = -2/(pi)*sin(2*pi*x);
sinf3 = 2/(3*pi)*sin(3*pi*x);
sinf = sinf1+sinf2+sinf3;

cosf1 = -2/pi*cos(pi*x);
cosf2 = 2/(3*pi)*cos(3*pi*x);
cosf3 = -2/(5*pi)*cos(5*pi*x);
cosf = 1/2+cosf1+cosf2+cosf3;

intf = (8/pi^3)*(sin(pi*x)+1/27*sin(3*pi*x))+2/pi*(sin(pi*x)-1/2*sin(2*pi*x));

figure(1)
plot(x,cosf)
title('Cosine Fourier Series')
xlabel('x')
ylabel('C(x)')

figure(2)
plot(x,sinf)
title('Sine Fourier Series')
xlabel('x')
ylabel('S(x)')

figure(3)
plot(x,intf)
title('Integration Fourier Series')
xlabel('x')
ylabel('S(x)')

