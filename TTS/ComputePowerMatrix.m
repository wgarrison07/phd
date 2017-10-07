function powerMatrix = ComputePowerMatrix(Order, InputDim)
    %Produces a Matrix of Size nxInputDim, where each row is a polynomial
    %variable and each column is the power to which the input variable is
    %raised.
    powerMatrix = zeros((Order+1)^InputDim, InputDim);
    for i = 2:size(powerMatrix,1)
        powerMatrix(i,:) = powerMatrix(i-1,:);
        for j = 1:InputDim
            if powerMatrix(i,j) == Order
                powerMatrix(i,j) = 0;
            else
                powerMatrix(i,j) = powerMatrix(i,j) + 1;
                break;
            end
        end
    end
end