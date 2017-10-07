function value = Michalewicz(input)
    sum = 0;
    for i = 1:length(input)
       sum = sum + sin(input(i))*sin(i*input(i)^2/pi)^20;        
    end
    value = -sum;
end