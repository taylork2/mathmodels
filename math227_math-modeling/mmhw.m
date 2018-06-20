P1 = [[5/18 5/18 2/9]; [0 5/12 0]; [0 0 5/9]];
ones = [1 1 1];
I=diag(ones);
N=inv(I-P1)
t=zeros(3,1);
for i=1:3
    t(i,1)=N(i,1)+N(i,2)+N(i,3);
end
t
R=[[0 0 2/9 0]; [5/12 0 1/12 1/12];[0 5/18 1/9 1/18]];
B=N*R

%{
P1=zeros(7,7);

for i=1:6
    P1(i,i+1)=0.4;
    P1(i+1, i)=0.6;
end
ones = [1 1 1 1 1 1 1];
I=diag(ones);
N=inv(I-P1)
R=zeros(7,2);
R(1,2)=0.6;
R(7,1)=0.4;
R
B=N*R

P2=[[0 0 0.4]; [0 0 0]; [0 0.6 0]];
ones = [1 1 1];
I=diag(ones);
N=inv(I-P2);
R=[[0.6 0]; [0.6 0.4]; [0 0.4]];
B=N*R

P3=[[0 0 0.4]; [.6 0 .4]; [0 0 0]];
ones = [1 1 1];
I=diag(ones);
N=inv(I-P3);
R=[[0.6 0]; [0 0]; [0.6 0.4]];
B=N*R

P4=[[0 .4 0]; [.6 0 .4]; [0 0.6 0]];
ones = [1 1 1];
I=diag(ones);
N=inv(I-P4)
R=[[0.6 0]; [0 0]; [0 0.4]];
B=N*R

P5 = [[.25 0 .25 .25]; [0 0 0 1]; [1 0 0 0]; [.25 .25 0 .25]];
ones = [1 1 1 1];
I=diag(ones);
N=inv(I-P5)
R=[[.25 0]; [0 0]; [0 0];[0 .25]];
B=N*R
%}

