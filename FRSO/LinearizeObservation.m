function H = LinearizeObservation(pos, coef, powerMat)

stateSize = length(coef);
H = zeros(length(pos)+1, stateSize);
for i = 1:stateSize
    step = coef;
    step(i) = step(i) + 1E-14*1i;
    H(:,i) = imag(ComputeEstimate(pos, step, powerMat))./1E-14;
end    

end