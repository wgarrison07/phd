function val = Shekel(input)
persistent m B C
if isempty(m)
    m = 10;
    B = 0.1*[1, 2, 2, 4, 4, 6, 3, 7, 5, 5];
    C = [4, 1, 8, 6, 3, 2, 5, 8, 6, 7;
         4, 1, 8, 6, 7, 9, 3, 1, 2, 3;
         4, 1, 8, 6, 3, 2, 5, 8, 6, 7;
         4, 1, 8, 6, 7, 9, 3, 1, 2, 1];
end

val = 0;
for i = 1:m
    subval = 0;
    for j = 1:4
        subval = subval + (input(j) - C(j,i))^2 + B(i);
    end
    val = val - 1/subval;
end

end