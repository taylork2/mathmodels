clear;
close all;

L = 2*pi/22; mu = 1; rho = 1; k = 22; eps = .001; g = 0; theta = pi/3; n =3; m =2;
sigma = 1;

h0 = .01; h_s = h0/10;
dt = .01; t = 0:dt:2;
dx = .001; x = 0:dx:L;

s1 = rho/(3*mu); r1 = 1/(8*dx^4);
s2 = 1/(3*mu); r2 = 1/(8*dx^2);

h = h0*ones(length(x),1) + eps*h0*sin((k)*x)';
dPhi = @(h1) (rho*(1-cos(theta))/((1/2)*h_s))*(-3*h_s^3/(h1^(3+1))+2*h_s^2/(h1^(2+1)));


fPm2 = @(hm2,hm1,h,hp1,hp2) r1*s1*(h + hm1)^3;
fPm1 = @(hm2,hm1,h,hp1,hp2) r2*s2*((rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((h+hm1)/2)^(3+1))+2*h_s^2/(((h+hm1)/2)^(2+1)))*(h + hm1)^3 - 3*(rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((h+hm1)/2)^(3+1))+2*h_s^2/(((h+hm1)/2)^(2+1)))*(h + hm1)^2*(h - hm1)) - r1*s1*(3*(h + hm1)^3 - 3*(h + hm1)^2*(3*h - 3*hm1 + hm2 - hp1) + (h + hp1)^3);
fP = @(hm2,hm1,h,hp1,hp2) r1*s1*(3*(h + hm1)^2*(3*h - 3*hm1 + hm2 - hp1) + 3*(h + hp1)^2*(3*h - hm1 - 3*hp1 + hp2) + 3*(h + hm1)^3 + 3*(h + hp1)^3) - r2*s2*((rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((h+hm1)/2)^(3+1))+2*h_s^2/(((h+hm1)/2)^(2+1)))*(h + hm1)^3 + (rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((hp1+h)/2)^(3+1))+2*h_s^2/(((hp1+h)/2)^(2+1)))*(h + hp1)^3 + 3*(rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((h+hm1)/2)^(3+1))+2*h_s^2/(((h+hm1)/2)^(2+1)))*(h + hm1)^2*(h - hm1) + 3*(rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((hp1+h)/2)^(3+1))+2*h_s^2/(((hp1+h)/2)^(2+1)))*(h + hp1)^2*(h - hp1)) + 1;
fPp1 = @(hm2,hm1,h,hp1,hp2) r2*s2*((rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((hp1+h)/2)^(3+1))+2*h_s^2/(((hp1+h)/2)^(2+1)))*(h + hp1)^3 - 3*(rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((hp1+h)/2)^(3+1))+2*h_s^2/(((hp1+h)/2)^(2+1)))*(h + hp1)^2*(h - hp1)) - r1*s1*((h + hm1)^3 - 3*(h + hp1)^2*(3*h - hm1 - 3*hp1 + hp2) + 3*(h + hp1)^3);
fPp2 = @(hm2,hm1,h,hp1,hp2) r1*s1*(h + hp1)^3;

g = @(hm2,hm1,h,hp1,hp2) h + r1*s1*((h+hp1)^3*(hp2-3*hp1+3*h-hm1)-(h+hm1)^3*(hp1-3*h+3*hm1-hm2))+ s2*r2*((h+hp1)^3*(rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((hp1+h)/2)^(3+1))+2*h_s^2/(((hp1+h)/2)^(2+1)))*(hp1-h)-(h+hm1)^3*(rho*(1-cos(theta))/(mu*h_s))*(-3*h_s^3/(((h+hm1)/2)^(3+1))+2*h_s^2/(((h+hm1)/2)^(2+1)))*(h-hm1));
indexHeight = find(h == max(h),1);

H = zeros(1,length(t));
Y = zeros(1,length(t));


figure(1)
xlabel('x');
ylabel('height');
title('Implicit Height of water w/ van der waals');

for j = 1:length(t)
    
    disp(j);
    H(j)= h(indexHeight);
    Y(j) = real(log((H(j)-h0)/(eps*h0*10*exp(1i*k*(1/4)))));
    
%     if (1 == mod(j,5))
        figure(1)
        hold on;
        plot(x,h');
%         ylim([.0997 .1003]);

%     end
    
    x0 = h;
    
    for i = 1:50
        J = zeros(length(h), length(h));
        
        %constructing the jacobian matrix
        for l = 3:length(h)-2
            J(l,l-2)= fPm2(x0(l-2),x0(l-1), x0(l), x0(l+1),x0(l+2));
            J(l,l-1)= fPm1(x0(l-2),x0(l-1), x0(l), x0(l+1),x0(l+2));
            J(l,l)= fP(x0(l-2),x0(l-1), x0(l), x0(l+1),x0(l+2));
            J(l,l+1) = fPp1(x0(l-2),x0(l-1), x0(l), x0(l+1),x0(l+2));
            J(l,l+2) = fPp2(x0(l-2),x0(l-1), x0(l), x0(l+1),x0(l+2));
        end

        %fixing the first column and the last column of jacobian manually
        
        J(1,end-2) = fPm2(x0(end-2), x0(end-1), x0(1),x0(2), x0(3));
        J(1,end-1) = fPm1(x0(end-2), x0(end-1), x0(1),x0(2), x0(3));
        J(1,1) = fP(x0(end-2), x0(end-1), x0(1),x0(2), x0(3));
        J(1,2)= fPp1(x0(end-2), x0(end-1), x0(1),x0(2), x0(3));
        J(1,3)= fPp2(x0(end-2), x0(end-1), x0(1),x0(2), x0(3));
        
        J(2,end-1) = fPm2(x0(end-1),x0(1),x0(2),x0(3),x0(4));
        J(2,1) = fPm1(x0(end-1),x0(1),x0(2),x0(3),x0(4));
        J(2,2) = fP(x0(end-1),x0(1),x0(2),x0(3),x0(4));
        J(2,3) = fPp1(x0(end-1),x0(1),x0(2),x0(3),x0(4));
        J(2,4) = fPp2(x0(end-1),x0(1),x0(2),x0(3),x0(4));

        J(end-1,2) = fPp2(x0(end-3),x0(end-2),x0(end-1),x0(end),x0(2));
        J(end-1,end) = fPp1(x0(end-3),x0(end-2),x0(end-1),x0(end),x0(2));
        J(end-1,end-1) = fP(x0(end-3),x0(end-2),x0(end-1),x0(end),x0(2));
        J(end-1,end-2) = fPm1(x0(end-3),x0(end-2),x0(end-1),x0(end),x0(2));
        J(end-1,end-3) = fPm2(x0(end-3),x0(end-2),x0(end-1),x0(end),x0(2));
        
        J(end,3) = fPp2(x0(end-2),x0(end-1),x0(end),x0(2),x0(3));
        J(end,2) = fPp1(x0(end-2),x0(end-1),x0(end),x0(2),x0(3));
        J(end,end) = fP(x0(end-2),x0(end-1),x0(end),x0(2),x0(3));
        J(end,end-1) = fPm1(x0(end-2),x0(end-1),x0(end),x0(2),x0(3));
        J(end,end-2) = fPm2(x0(end-2),x0(end-1),x0(end),x0(2),x0(3));
        
        J = sparse(J);
        
        
        F = zeros(length(h),1);
        
        for k = 3:length(h)-2
            % contructing f(x) from newton's method X_n+1 = X_n + [J]^-1*f(x)
            F(k) = -h(k) + g(x0(k-2),x0(k-1), x0(k), x0(k+1),x0(k+2));
        end
        
        %fixing the first point and the last point
        F(1) = -h(1) + g(x0(end-2),x0(end-1), x0(1), x0(2),x0(3));
        F(2) = -h(2) + g(x0(end-1),x0(1), x0(2), x0(3),x0(4));
        F(end-1) = -h(end-1) + g(x0(end-3),x0(end-2),x0(end-1),x0(end),x0(2));
        F(end) = -h(end) + g(x0(end-2),x0(end-1),x0(end),x0(2),x0(3));
        
        h_new = x0 - J\F;
        
        if (norm(x0 - h_new) < .0000001) %use the 2 norm for the error
            break; 
        end
        
        
        x0 = h_new;
    end
    
    h = x0;
end


w = @(k) -(h0^3/(3*(3-2)/((2-1)*(3-1))))*(sigma.*k.^4 - k.^2.*sigma.*(1-.7071)./((1/2).*h_s).*((-3.*h_s.^3)./h0.^4+ (2.*h_s.^2)./h0.^3));

figure 
plot(x, w(x));
figure
title('Slope at x_{max} vs t');
plot(Y(1:10));
xlabel("t");
ylabel("y");
