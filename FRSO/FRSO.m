function [optInput, optVal, vectorField] = FRSO(objFunc, varRanges, varargin)
%% Input Error Checking
if ~isa(objFunc, 'function_handle')
    error('objFunc must be a function handle');
end

if size(varRanges,1) ~= 2
    error('varRanges must be 2xn in size');
end


%% Initialize Variables
inputDim = size(varRanges,2);
maxOrder = max(6 - inputDim, 2);
maxIter = 5000;
modelErr = zeros(maxIter, 1);
iIter = 0;
reqdPopt = 0.95;
debugPlot = true;
priming = true;
numPrimingPoints = ceil(100*sqrt(inputDim));
primingPoints = RandomLatinHypercube(varRanges, numPrimingPoints);
strLen = 0;
sumErr = 0;
sumVal = 0;
sumWgt = 0;

%% Process Input Arguments
iArg = 1;
nArg = length(varargin);
while iArg <= nArg
    iArg = iArg + 1;
end

%% Create Model Coeffcients
numCoef = (maxOrder + 1)^inputDim;
coef = zeros(numCoef, 1);
coefCov = 1E10*eye(numCoef);
powerMat = Permutations(numCoef, inputDim, maxOrder);

%% Main Loop
nextPos = primingPoints(iIter+1,:)';
while iIter < maxIter
    %Update Iteration Counter and Print Status to Screen
    iIter = iIter + 1;
    fprintf(1, repmat('\b',1,strLen));
    strLen = fprintf(1, 'Executing Iteration: %i', iIter);
        
    %Evaluate next Position
    [valVect, valCov] = EvalObjFunc(objFunc, nextPos, 1);
    
    %Filter Response Surface
    y = valVect - ComputeEstimate(nextPos, coef, powerMat);
    H = LinearizeObservation(nextPos, coef, powerMat);
    S = H*coefCov*H' + valCov;
    K = coefCov*H'/S;
    coef = coef + K*y;
    ImKH = eye(length(coef)) - K*H;
    coefCov = ImKH*coefCov*ImKH' + K*valCov*K';
    inflator = 1.05*eye(size(coefCov));
    coefCov = inflator*coefCov*inflator';
    
    %Compute Model Accuracy Statistics
    sumErr = sumErr + iIter^2*abs(y(1));
    sumVal = sumVal + abs(valVect(1));
    sumWgt = sumWgt + iIter^2;
    modelErr(iIter) = (sumErr/sumWgt)/(sumVal/iIter);
    
    %Debug Plotting if Needed
    if debugPlot && iIter > 220
        PlotSurf(varRanges, coef, powerMat);
        bob = 1;
    end
    
    %Compute Next Position
    if priming
        nextPos = primingPoints(iIter+1,:);
        if iIter+1 == numPrimingPoints
            priming = false;
            coef = [coef; zeros(size(coef))];
            newCov = 1E8*eye(length(coef));
            newCov(1:numCoef, 1:numCoef) = coefCov;
            coefCov = newCov;
            powerMat = [powerMat; (powerMat + maxOrder)];
        end
    else
        nextPos = [0, 0] + randn(1,2)*0.1;
        test = 2*eye(size(coefCov));
        coefCov = test*coefCov*test';
        %[nextPos, score, primingPoints] = SelectNextPosition(varRanges, primingPoints, coef, powerMat, coefCov);
    end
    
end

end


function PlotSurf(varRanges, coef, powerMat)
    x = varRanges(1,1):0.1:varRanges(2,1);
    y = varRanges(1,2):0.1:varRanges(2,2);
%     x = -0.5:0.05:0.5;
%     y = -0.5:0.05:0.5;
    [X, Y] = meshgrid(x, y);
%     Z = X.^2 + Y.^2;
%     Z = (1 - X.^2) + 100*(Y - X.^2).^2;
    Z = -20*exp(-0.2*sqrt(0.5*(X.^2 + Y.^2))) - exp(0.5*(cos(2*pi()*X) + cos(2*pi()*Y))) + 20;
%     Z = 2*X.^2 - 1.05*X.^4 + X.^6/6 + X.*Y + Y.^2;
    
    gcf();
    hold off;
    mesh(X, Y, Z);
    hold on;
    for i = 1:length(x)
        for j = 1:length(y)
            val = ComputeEstimate([x(i), y(j)], coef, powerMat);
            plot3(x(i), y(j), val(1) , '.g');
        end
    end
    xlabel('x1');
    ylabel('x2');
end


function permMat = Permutations(n, k, max)
permMat = zeros(n, k);
line = 1;
perm = zeros(1, k);
while line <= n
    permMat(line,:) = perm;
    for i = 1:k
        if perm(i) == max
            perm(i) = 0;
        else
            perm(i) = perm(i) + 1;
            break;
        end
    end
    line = line + 1;
end
end