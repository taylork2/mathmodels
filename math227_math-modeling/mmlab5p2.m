N=100000;
h8=0;
h20=0;
sumHours = 0;
for i=1:N
    state = 1;
    hours = 0;
    while (state ~= 3)
        r=rand;
        hours=hours+1;
        if (state == 1)
            if (r < 1/3)
               state = 1;
            elseif (r > 2/3)
                state = 2;
            else 
                state = 3;
            end
        elseif (state == 2)
            if (r <1/2)
                state = 1;
            else
                state = 2;
            end
        end
    end
    sumHours = sumHours + hours;
    if (hours <= 8) 
        h8 = h8+1;
    elseif (hours > 20)
        h20 = h20+1;
    end
end

hrsEscape = sumHours/N
prob8 = h8/N
prob20 = h20/N
            
        