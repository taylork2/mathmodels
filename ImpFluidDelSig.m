clear;
%declare constants
p=1;  g=0; mu=1;  eps = 0.1; h0 = 1; k = 2*pi; sig = 1; delsig = 0; alpha_th = 1; k_th = 1;
%declare spacing in the x-direction
x = linspace(0,1,60);    dx = x(3)-x(2);    N = length(x);
%declare time step spacing
t = linspace(0,1,100);            dt = t(2)-t(1);
%Make tolerance smaller, for more iterations of Newton's method
tol = 0.1;
  
S = (dt)/(24*mu*(dx^2));
a1 = p*g;       %constant infront of gravity scheme
a2 = sig/(dx^2); %contant in front of surface tension scheme
a3 = alpha_th*delsig/(2*k_th)*(1/(4*dx^2));

%declare initial conditions and plot
h=zeros(N,length(t));   %h_max=zeros(length(t));
for i = 1:1:N
    h(i,1)=h0+(h0*eps*sin(k*x(i)) );
end
hmax(1) = max(h(:,1));
h_max(1) = h(floor(N/4),1);
oldh = h;    newh = h;
    subplot(2,1,1); plot(x,h(:,1),'k-','Linewidth',1.5),    hold on;
    xlabel('X'); ylabel('H'); title('Fluid Equation with Marangoni Effect');
    ylim([h0-1.1*eps,h0+1.1*eps]); 

%Implicit scheme as a function  
%h(i-2) = x1    h(i-1) = x2     h(i) = x3   h(i+1)=x4  h(i+2) = x5
%Gravity scheme (original)
g1 = @(x2,x3,x4) ((x4+x3)^3) *(x4-x3) - ((x3+x2)^3)*(x3-x2);
%Surface Tension scheme
g2 = @(x1,x2,x3,x4,x5) ((x4+x3)^3) * (x5-3*x4+3*x3-x2) - ((x3+x2)^3) * (x4-3*x3+3*x2-x1);

g3 = @(x1,x2,x3,x4,x5) (x4-x3) * (x4+x3)^2 - (x3-x2)*(x3+x2)^2;

%Total scheme with Gravity and Surface Tension
f= @(x1,x2,x3,x4,x5,oldx) x3-S*( a1*g1(x2,x3,x4) - a2*g2(x1,x2,x3,x4,x5) - a3*g3(x1,x2,x3,x4,x5)) - oldx;       

%f= @(x1,x2,x3,x4,x5,oldx) x3-S*( a1*((x4+x3)^3*(x4-x3)-((x3+x2)^3*(x3-x2))) - a2*((x4+x3)^3*(x5-3*x4+3*x3-x2)-((x3+x2)^3*(x4-3*x3+3*x2-x1)) ) - a3*((x4-x3) * (x4+x3)^2 - (x3-x2)*(x3+x2)^2) ) -oldx;       

%derivatives of the scheme to use in the Jacobian
%derivative wrt x1
df1=@(x1,x2,x3,x4,x5) S*a2*(x2 + x3)^3;  
%derivative wrt x2
df2=@(x1,x2,x3,x4,x5) -S*(a2*(3*(x2 + x3)^3 - 3*(x2 + x3)^2*(x1 - 3*x2 + 3*x3 - x4) + (x3 + x4)^3) + a1*((x2 + x3)^3 + 3*(x2 + x3)^2*(x2 - x3)) - a3*((2*x2 + 2*x3)*(x2 - x3) + (x2 + x3)^2)); 
%derivative wrt x3    
df3=@(x1,x2,x3,x4,x5) S*(a2*(3*(x2 + x3)^2*(x1 - 3*x2 + 3*x3 - x4) - 3*(x3 + x4)^2*(x2 - 3*x3 + 3*x4 - x5) + 3*(x2 + x3)^3 + 3*(x3 + x4)^3) - a3*((2*x3 + 2*x4)*(x3 - x4) - (2*x2 + 2*x3)*(x2 - x3) + (x2 + x3)^2 + (x3 + x4)^2) + a1*((x2 + x3)^3 + (x3 + x4)^3 - 3*(x2 + x3)^2*(x2 - x3) + 3*(x3 + x4)^2*(x3 - x4))) + 1;  
%derivative wrt x4    
df4=@(x1,x2,x3,x4,x5) -S*(a2*(3*(x3 + x4)^2*(x2 - 3*x3 + 3*x4 - x5) + (x2 + x3)^3 + 3*(x3 + x4)^3) + a3*((2*x3 + 2*x4)*(x3 - x4) - (x3 + x4)^2) + a1*((x3 + x4)^3 - 3*(x3 + x4)^2*(x3 - x4)));
%derivative wrt x5    
df5=@(x1,x2,x3,x4,x5) S*a2*(x3 + x4)^3;

J = zeros(N-2,N-2);
b=zeros(N-2,1);     

c=1;    %variable used for plotting
for tstep = 1:length(t)-1
      %make the jacobian    
    for j = 0:tol:1 
            %First row
            J(1,1)    = df3(newh(N-1,tstep),newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep));
            J(1,2)    = df4(newh(N-1,tstep),newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep)); 
            J(1,3)    = df5(newh(N-1,tstep),newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep)); 
            %second row
            J(2,1)    = df2(newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep),newh(5,tstep));
            J(2,2)    = df3(newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep),newh(5,tstep)); 
            J(2,3)    = df4(newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep),newh(5,tstep)); 
            J(2,4)    = df5(newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep),newh(5,tstep)); 
            %second to last row
            J(N-3,N-5)      = df1(newh(N-4,tstep),newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep));   
            J(N-3,N-4)      = df2(newh(N-4,tstep),newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep));
            J(N-3,N-3)      = df3(newh(N-4,tstep),newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep));  
            J(N-3,N-2)      = df4(newh(N-4,tstep),newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep));
            %last row
            J(N-2,N-4)      = df1(newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep),newh(1,tstep));   
            J(N-2,N-3)      = df2(newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep),newh(1,tstep));
            J(N-2,N-2)      = df3(newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep),newh(1,tstep));
            
        for i = 3:N-4   %Middle Rows
            J(i,i-2)    = df1(newh(i-1,tstep),newh(i,tstep),newh(i+1,tstep),newh(i+2,tstep),newh(i+3,tstep));   
            J(i,i-1)    = df2(newh(i-1,tstep),newh(i,tstep),newh(i+1,tstep),newh(i+2,tstep),newh(i+3,tstep));            
            J(i,i)      = df3(newh(i-1,tstep),newh(i,tstep),newh(i+1,tstep),newh(i+2,tstep),newh(i+3,tstep));
            J(i,i+1)    = df4(newh(i-1,tstep),newh(i,tstep),newh(i+1,tstep),newh(i+2,tstep),newh(i+3,tstep)); 
            J(i,i+2)    = df5(newh(i-1,tstep),newh(i,tstep),newh(i+1,tstep),newh(i+2,tstep),newh(i+3,tstep)); 
        end 
        
         %make the b vector
        b(1,tstep) = f(newh(N-1,tstep),newh(1,tstep),newh(2,tstep),newh(3,tstep),newh(4,tstep),oldh(2,tstep));
        b(N-2,tstep) = f(newh(N-3,tstep),newh(N-2,tstep),newh(N-1,tstep),newh(N,tstep),newh(2,tstep),oldh(N-1,tstep)); 
        for i = 2:N-3     
             b(i,tstep) = f(newh(i-1,tstep),newh(i,tstep),newh(i+1,tstep),newh(i+2,tstep),newh(i+3,tstep),oldh(i+1,tstep));
        end   
        
        %computation of newton's method
         backslash = J\b(:,tstep);       %seperated math for troubleshooting 
         newh(2:N-1,tstep) = newh(2:N-1,tstep) - backslash;
         newh(1,tstep) = h0; newh(N,tstep) = h0;
    end %newtons method loop end
    %update the oldh and h
      oldh(2:N-1, tstep+1) = h(2:N-1,tstep);    
      h(1:N,tstep+1) = newh(1:N,tstep);
       hmax(tstep+1) = max(h(:,tstep+1));
       h_max(tstep+1) = h(floor(N/4),tstep+1);
     %plotting iterations
        if c == 1
        subplot(2,1,1); plot(x,newh(:,tstep)); hold on;
        end
        if c >=5   %change this value to a higher number to plot less iterations
            c = 0;
        end
        c=c+1; %updates c every iteration
    
end

%plot the last iteration in case it was not plotted in the loop
subplot(2,1,1)
plot(x,h(:,length(t)),'k-','Linewidth',1.5);

%Check the slope using max h values
    F = log(abs(h_max - h0)) - log(eps*h0);
    w = (k^2/mu*(alpha_th*delsig*(h0)/(k_th*2) -(h0).^3*(sig*k^2 + p*g)/3))/length(t)
   
    subplot(2,1,2)
    plot(0:dt:1, F);    % xlim([0,50])       
    xlabel('t'); ylabel('w(t)'); title('Slope Checker');
    
    final = floor(length(t)/2);
    fluidslope = (h_max(final) - h_max(1))/ (t(final) - t(1));    
    disp(['Fluid Experiemental Rate of Decay: ', num2str(fluidslope)]);