function value = Ackley(input)

a = 20;
b = 0.2;
c = 2*pi;
invDim = 1/length(input);
value = -a*exp(-b*sqrt(sum(input.^2)*invDim)) - exp(sum(cos(c*input))*invDim) + a + exp(1);

end