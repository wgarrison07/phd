function value = Eggholder(input)

x = input(1);
y = input(2);

value = -(y+47)*sin(sqrt(abs(x/2+y+47))) - x*sin(sqrt(abs(x - y - 47)));

end