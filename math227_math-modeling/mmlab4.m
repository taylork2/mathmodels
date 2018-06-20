S=table3;
S=flip(S);
days = size(S,1);
h=1/days;
X=1:(days-1);

for i = 1:days-1
     X(i)= log(S(i+1)/S(i));
end

v = mean(X)
T = var(X)  

sigma = sqrt(T)/sqrt(h) %adjusted standard deviation
u = v./h + T./(2*h) %adjusted mean

N=1000;
Sto = 1:days;
Sto(1) = S(days);
m=1:N;
for n=1:N
    index=1;
    for t=0:days
        Sto(index+1) = Sto(index)*exp(h*(u-sigma^2/2)+sigma*normrnd(0, sqrt(h)));
        index=index+1;
    end  
    
    m(n)=Sto(days);
end
hist(m)
prob = 1-normcdf(log(1.1)/days, u, sigma)
prob2 = 1-logncdf(1.1, u, sigma)
