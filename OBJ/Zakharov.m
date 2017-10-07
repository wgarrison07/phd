function value = Zakharov(input)

sum1 = 0;
sum2 = 0;
for i = 1:length(input)
    sum1 = sum1 + input(i)^2;
    sum2 = sum2 + 0.5*i*input(i);
end

value = sum1 + sum2^2 + sum2^4;

end