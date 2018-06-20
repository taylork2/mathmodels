N=100;
M = zeros(N,2);

for i=1:N-1
    r=rand;
    if (r<0.25)
        M(i+1,1)=M(i,1);
        M(i+1,2)=M(i,2)+1;
    elseif (r<0.5)
        M(i+1,1)=M(i,1)+1;
        M(i+1,2)=M(i,2);
    elseif (r <0.75)
        M(i+1,1)=M(i,1);
        M(i+1,2)=M(i,2)-1;
    else
        M(i+1,1)=M(i,1)-1;
        M(i+1,2)=M(i,2);
    end 
end 
plot(M(:,1), M(:,2))
title 'Positions at each step'

away=0;
k=1000;
for j=1:k
    for i=1:N-1
        r=rand;
        if (r<0.25)
            M(i+1,1)=M(i,1);
            M(i+1,2)=M(i,2)+1;
        elseif (r<0.5)
            M(i+1,1)=M(i,1)+1;
        	M(i+1,2)=M(i,2);
        elseif (r <0.75)
            M(i+1,1)=M(i,1);
            M(i+1,2)=M(i,2)-1;
        else
            M(i+1,1)=M(i,1)-1;
            M(i+1,2)=M(i,2);
        end
    end 
    if abs(M(i+1,1))>=10 || abs(M(i+1,2))>=10
        away=away+1;
    end
end 
prob=away/k
            