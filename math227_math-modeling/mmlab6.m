%{
array = size(mmlab6p1);
N = array(1,1);
hold on
for i=1:N
    plot(mmlab6p1(i, 1), mmlab6p1(i,2), '*')
end
hold off
%}


plot(log(BodyWeightg), log(Pulseratebeatsmin))
xlabel('Body Weight (g)')
ylabel('Pulse rate (beats/min)')
title ('Pulse rate vs Body Weight')

length = size(VarName2);
years = 1:length(1,1);
plot(years, log(VarName2))
xlabel('Year')
ylabel('Population')
title('Population Over Years')
 