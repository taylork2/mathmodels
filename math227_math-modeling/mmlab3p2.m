clear all;
h=0.01;
u=0.1;
o=0.5;
S0=1.0;
S=0:0.01:5;
S(1)=S0;
m=1:10000;

for n=1:10000
    index=1;
    for t=0:0.01:5
        S(index+1) = S(index)*exp(h*(u-o^2/2)+o*normrnd(0, sqrt(h)));
        index=index+1;
    end  
    
    m(n)=S(502);
end
Mean = mean(m)
ConfInt = [Mean-1.96*o/sqrt(10000), Mean+1.96*o/sqrt(10000)]
hist(log(m))
title('Histogram of log(S_{500,i})')
xlabel('dS_t/S_t')
ylabel('Frequency')