

function model = ANN(x, y, varargin)

%% Input Processing
if nargin < 2
    error('x and y are required inputs');
elseif size(x,1) ~= size(y,1)
    error('Fitting Inputs and Outputs must have the same size');
end

%% Initialize Values
model = SurrogateModel();
model.Type = 'ANN';
model.InputDim = size(x,2);
model.OutputDim = size(y,2);
numSamples = size(x,1);
maxStallCount = 100;
testPlot = true;

%Normalize and Sort Inputs, Normalize Outputs
[x, model.InputMu, model.InputSigma] = NormalizeData(x);
[y, model.OutputMu, model.OutputSigma] = NormalizeData(y);
sumX = sum(x,2);
[~, indx] = sort(sumX);
x = x(indx,:);
y = y(indx,:);

%% Divide Data into Training and Testing Sets
testP     = 0.2; %Twenty Percent
trainIndx = mod(1:numSamples, 1/testP) ~= 0;
trainX    = x(trainIndx, :);
trainY    = y(trainIndx, :);
testX     = x(~trainIndx, :);
testY     = y(~trainIndx, :);

%% Initialize Neural Network
model.numHLN = 10 + 2*ceil(sqrt(size(trainX,1)));
nodes = (1:model.numHLN)';
W1 = randn(model.numHLN, model.InputDim + 1);
W1(:,end-1) = (nodes - mean(nodes))/std(nodes); 
W2 = randn(model.OutputDim, model.numHLN + 1);
model.MVE = Inf;
stallCount = 0;

%% Model Training Loop
alpha = 1;
for iter = 0:1E4
    %Compute Backpropagation Error Gradients
    [G1, G2] = ComputeErrorGradient(W1, W2, trainX, trainY);
    
    %Correct Neuron Weights
    W1 = W1 - G1*alpha;% - H1*alpha^2;
    W2 = W2 - G2*alpha;% - H2*alpha^2;
    
    %Compute Testing Error
    yEst = Predict(W1, W2, testX);
    MVE = mean(sqrt(sum((yEst - testY).^2, 2)));
    if mod(iter, 100) == 0 || true
        fprintf(1, 'iter = %f, MVE = %f, alpha = %f\n', iter, MVE, alpha);
    end
    
    %Check for Training Abort
    if MVE >= model.MVE
        alpha = max(0.9*alpha, 1E-9);
        W1 = model.W1;
        W2 = model.W2;
        stallCount = stallCount + 1;
        if stallCount > maxStallCount
            break;
        end
    else
       %Update Model
       alpha = alpha*(1 + 1E-4);
       stallCount = 0;
       model.MVE = MVE;
       model.W1 = W1;
       model.W2 = W2;
    end
end

%% Compute Fit Statistics

%Test Plot if Testing Algorithm
if (testPlot && true)
    if model.InputDim == 1
        figure();
        plot(x, y, '.k');
        hold on;
        plot(testX, yEst, '.g');
        yEst = Predict(W1, W2, trainX);
        plot(trainX, yEst, '.g');
    elseif model.InputDim == 2
        figure();
        plot3(x(:,1), x(:,2), y, '.k');
        hold on;
        plot3(testX(:,1), testX(:,2), yEst, '.g');
        yEst = Predict(W1, W2, trainX);
        plot3(trainX(:,1), trainX(:,2), yEst, '.g');
    end
end

end

function yEst = Predict(W1, W2, x)
    n = size(x, 1);
    yEst = zeros(n, size(W2, 1));
    x = [x, ones(n,1)];
    for i = 1:n
        z = W1*x(i, :)';
        g = [tanh(z);1];
        yEst(i, :) = (W2*g)';
    end

end

function [G1, G2] = ComputeErrorGradient(W1, W2, x, y) 

G1 = zeros(size(W1));
G2 = zeros(size(W2));
n = size(x, 1);
m = size(x, 2) + 1;
x = [x, ones(n,1)];
for i = 1:n
    %Forward Pass
    z = W1*x(i, :)';
    g = [tanh(z);1];
    yEst = (W2*g)';
    
    %Backward Pass
%     err = ((-1-y(i,:))./(1+yEst) - (-1+y(i,:))./(1-yEst));
    err = 2*(yEst - y(i,:));
    G2 = G2 + (g*err)';
    dedg = repmat(err*W2, m, 1)';
    grad = dedg.*((1 - g.^2)*x(i, :));
    G1 = G1 + grad(1:end-1,:);
end
G1 = G1/n;
G2 = G2/n;
    
end
