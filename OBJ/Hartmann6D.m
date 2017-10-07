function val = Hartmann6D(input)
persistent alpha A P
if isempty(alpha)
    alpha = [1.0, 1.2, 3.0, 3.2];
    A = [10, 3, 17, 3.50, 1.7, 8;
         0.05, 10, 17, 0.1 8, 14;
         3, 3.5, 1.7, 10, 17, 8;
         17, 8, 0.05, 10, 0.1, 14];
     P = 1E-4*[1312, 1696, 5569, 0124, 8283, 5886;
               2329, 4135, 8307, 3736, 1004, 9991;
               2348, 1451, 3522, 2883, 3047, 6650;
               4047, 8828, 8732, 5743, 1091, 0381];
end
       
val = 0;
for i = 1:4
    subval = 0;
    for j = 1:6
        subval = subval + A(i, j)*(input(j) - P(i, j))^2;
    end
    val = val - alpha(i)*exp(subval);
end

end