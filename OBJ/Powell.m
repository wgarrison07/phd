function val = Powell(input)

val = 0;
for i = 1:length(input)/4
    val = val + (input(4*i-3) + 10*input(4*i-2))^2 + 5*(input(4*i-1) - input(4*i))^2 + (input(4*i-2) - 2*input(4*i-1))^4 + 10*(input(4*i-3) - input(4*i))^4;
end

end