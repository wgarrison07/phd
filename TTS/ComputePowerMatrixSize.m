function size = ComputePowerMatrixSize(Order, InputDim)

%Build Pascal's Triangle
coefs = zeros(Order);
coefs(:, 1) = 1;
for i = 2:Order
    for j = 2:Order
        coefs(i, j) = coefs(i-1, j) + coefs(i-1, j-1);
    end
end

size = 0;
for i = 1:Order
    if i > InputDim
        size = size + 1;
    else
        size = size + sum(coefs(:, i))*nchoosek(InputDim, i);
    end
end

end