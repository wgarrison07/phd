function [measEst, measEstCov] = PolyRegMeasFunc(PolyCoef, PolyCoefCov, X, PolyPowerMat)
%PolyCoef is a vector of size [c, 1] representing the polynomial coefficients
%PolyCoefCov is a [c, c] covariance matrix of the coefficients
%X is a [n, 1] vector of the desired input point
%PolyPowerMat is an [c, n] matrix of the power to which each dimesion is
%raised for each coefficient

    %Generate Vandermoode Matrix from Power Matrix
    numCoef = size(P,1);
    M = ones(1, numCoef);
    for i = 1:numCoef
        for j = 1:model.InputDim
            M(i) = M(i)*(X(j)^PolyPowerMat(i,j));
        end
    end
    
    %Compute estimated measurement
    measEst = M*PolyCoef;
    
    %Compute estimated measurement covariance
    measEstCov = M*PolyCoefCov*M';  

end