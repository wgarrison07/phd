function [nextPos, score, positions] = SelectNextPosition(varRange, positions, coef, powerMat, coefCov)
%Initialize Variables
numPoints = size(positions, 1);
fVals = zeros(numPoints,1);

%Perform Line Search using Quasi-Newton Search on All Points
for i = 1:numPoints
    estimate = ComputeEstimate(positions(i,:), coef, powerMat);
    s = -estimate(2:end)';
    s = s/norm(s);
    maxStep = ComputeStepLimit(positions(i,:), varRange, s);
    [step, fVals(i)] = fminbnd(@(x) LineSearchFunc(x, positions(i,:), s, coef, powerMat), 0, maxStep);
    positions(i, :) = positions(i, :) + step*s;
end

figure();
plot(sort(fVals));
end

function maxStep = ComputeStepLimit(pos, varRange, searchDir)
    maxStep = Inf;
    for i = 1:length(pos)
       if searchDir(i) > 0
          step = (varRange(2,i) - pos(i))/searchDir(i);
       else
          step = (varRange(1,i) - pos(i))/searchDir(i);
       end
       maxStep = min(step, maxStep);
    end
end

function val = LineSearchFunc(step, pos, searchDir, coef, powerMat)
    pos = pos + step*searchDir;
    valVect = ComputeEstimate(pos, coef, powerMat);
    val = valVect(1);
end

function [grad, hess] = PositionSensitivity(pos, coef, powerMat)

for i = 1:stateSize
    step = coef;
    step(i) = step(i) + 1E-14*1i;
    H(:,i) = imag(ComputeEstimate(pos, step, powerMat))./1E-14;
end    


end