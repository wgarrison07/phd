function value = Beale(input)
    x1 = input(1);
    x2 = input(2);
    value = (1.5 - x1 + x1*x2)^2 + (2.25 - x1 + x1*x2^2)^2 + (2.625 - x1 + x1*x2^3)^2;
end