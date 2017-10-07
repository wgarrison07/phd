function [optInput, optVal, vectorField] = FVF(objFunc, varRanges, varargin)
%% Input Error Checking
if ~isa(objFunc, 'function_handle')
    error('objFunc must be a function handle');
end

if size(varRanges,1) ~= 2
    error('varRanges must be 2xn in size');
end


%% Initialize Variables
inputDim = size(varRanges,2);
numNodes = ceil(sqrt(inputDim))+6;
eval2UpdateMesh = Inf;
maxIter = 5000;
iIter = 0;
reqdPopt = 0.95;
debugPlot = true;

%% Process Input Arguments
iArg = 1;
nArg = length(varargin);
while iArg <= nArg
    iArg = iArg + 1;
end

%% Initialize VectorField
vectorField = CreateVectorField(varRanges, numNodes);
nextPos = mean(varRanges,1);

%% Main Loop
while iIter < maxIter
    %Evaluate next Position
    iIter = iIter + 1;
    [valVect, valCov] = EvalObjFunc(objFunc, nextPos, 3);
    
    %Update Mesh Filter
    vectorField = FilterField(vectorField, nextPos, valVect, valCov);
    
    %Debug Plotting if Needed
    if debugPlot
        PlotFVF(varRanges, vectorField, nextPos, valVect);
    end
    
    %Compute P(optimum), Current Optimum, and Next Position to Evaluate
    [Popt, optInput, optVal, nextPos] = EvaluateField(vectorField);
    if Popt > reqdPopt
        break;
    end
    
    %Update Mesh Locations if Needed
    if mod(iIter, eval2UpdateMesh) == 0
        vectorField = UpdateField(vectorField);
    end
    
end


%% Set Outputs

end

function PlotFVF(varRanges, field, nextPos, valVect)
    x = varRanges(1,1):0.1:varRanges(2,1);
    y = varRanges(1,2):0.1:varRanges(2,2);
    [X, Y] = meshgrid(x, y);
    Z = X.^2 + Y.^2;
    gcf();
    hold off;
    mesh(X, Y, Z);
    hold on;
    plot3(field.nodePos(:,1), field.nodePos(:,2), field.nodeVal(:,1), '.k');
    for i = 1:length(x)
        for j = 1:length(y)
            val = ComputeEstimateFromField(field, [x(i), y(j)]);
            plot3(x(i), y(j), val(1) , '.g');
        end
    end

end


function [Popt, optInput, optVal, nextPos] = EvaluateField(field)
Popt = 0;
optInput = 0;
optVal = 0;
nextPos = rand(1,field.inputDim);
end