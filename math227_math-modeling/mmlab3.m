h=0.01;
u=0.1;
o=0.5;
S0=1.0;
S=0:0.01:5;
S(1)=S0;
n = 0:0.01:5.01;

for x=1:4
    index = 1;
    subplot(2,2,x);
    xlabel('Time (years)')
    ylabel('Stock Price ($)')
    title('Stock Price over Time')
    hold on
    for t=0:0.01:5
        S(index+1) = S(index)*exp(h*(u-o^2/2)+o*normrnd(0, sqrt(h)));
        index=index+1;
    end    
    plot(n,S)
end

hold off