x_array = 1:10;
for i = 1:10
    if mod(x_array(i), 3) == 1
        disp(x_array(i))
    elseif mod(x_array(i), 3) == 2 
        disp(-x_array(i) )
    else 
        disp(x_array(i)*10)
    end
end

