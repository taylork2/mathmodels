%{
Yexp = labn2(:,2); %distance converted to meter
Iexp = labn2(:,1); %light intensity
[v,ind] = max(Iexp) 
IexpRel = Iexp/max(Iexp);
YexpRe = (Yexp-Yexp(ind))*10^-2; 
plot(YexpRe,IexpRel)
xlabel('distance(m)')
ylabel('I/Imax')
axis([-0.05 0.05 0 1])
%The min and max value of the recalculated distance
min(YexpRe)
max(YexpRe)

L = .675;
a = 0.04*10^-3;
d= 0.25*10^-3;
lamda = 650*10^-9;
y = -0.05:.0001:0.05;
theta = atan(y./L);
alpha = pi*a*sin(theta)./lamda;
beta = pi *d*sin(theta)./lamda;
IratioDouble = (cos(beta).*sin(alpha)./alpha).^2;
plot(y, IratioDouble, ':', YexpRe, IexpRel)
xlabel 'distance(m)', ylabel 'I/Imax'
%}

L = .675;
a = 0.04*10^-3;
lamda = 100*10^-9;
PD = 2*pi*a*y./(lamda*L);
Iratiosingle = (sin(.5*PD)./(.5*PD)).^2;
plot(y,Iratiosingle,':')
xlabel 'distance(m)', ylabel 'I/Imax'
