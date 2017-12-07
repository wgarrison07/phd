function powerMatrix = ComputePowerMatrix(Order, InputDim)
    %Produces a Matrix of Size nxInputDim, where each row is a polynomial
    %variable and each column is the power to which the input variable is
    %raised.
    size = ComputePowerMatrixSize(Order, InputDim);
    powerMatrix = zeros(size, InputDim);
    for i = 2:size
        powerMatrix(i,:) = powerMatrix(i-1,:);
        for j = 1:InputDim
            if sum(powerMatrix(i,:)) == Order
                powerMatrix(i, j) = 0;
            else
                powerMatrix(i, j) = powerMatrix(i, j) + 1;
                break;
            end
        end
    end

end