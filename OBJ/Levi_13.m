function value = Levi_13(input)

x = input(1);
y = input(2);
value = sin(3*pi*x)^2 + (x-1)^2*(1+sin(3*pi*y)^2) + (y-1)^2*(1+sin(2*pi*y)^2);

end