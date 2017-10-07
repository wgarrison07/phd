function value = Bukin_6(input)
    x1 = input(1);
    x2 = input(2);
    value = 100*sqrt(abs(x2 - 0.01*x1^2)) + 0.01*abs(x1 + 10);
end