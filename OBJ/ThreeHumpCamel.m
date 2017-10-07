function value = ThreeHumpCamel(input)

x = input(1);
y = input(2);
value = 2*x^2 -1.05*x^4 + x^6/6 + x*y + y^2;

end