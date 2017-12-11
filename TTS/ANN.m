

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
maxStallCount = 5;
testPlot = false;

%Normalize and Sort Inputs, Normalize Outputs
[x, model.InputMu, model.InputSigma] = NormalizeData(x);
[y, model.OutputMu, model.OutputSigma] = NormalizeData(y);
sumX = sum(x,2);
[~, indx] = sort(sumX);
x = x(indx,:);
y = y(indx,:);

%% Initialize Neural Network
model.numHLN = ceil(sqrt(size(x, 1)))*model.InputDim + 10;
maxAlp = 1/model.numHLN;
W1 = randn(model.numHLN, model.InputDim+1);
W2 = randn(model.OutputDim, model.numHLN + 1);
model.MVE = 1E8;
stallCount = 0;


%% Model Training Loop
for iter = 0:1E3
    %Compute Backpropagation Error Gradients
    [G1, G2] = ComputeErrorGradient(W1, W2, x, y);
    
    %Correct Hidden Layer Neuron Weights w/ LineSearch
    alpha = BoundLineSearch(@(alpha) Error(W1, W2, G1, G2, x, y, alpha), 0, maxAlp);
    W1 = W1 - alpha*G1;
    W2 = W2 - alpha*G2;
    
    %Compute Estimate Validation Error (training error)
    yEst = PredictANN(W1, W2, x);
    MVE = sqrt(mean(sum((yEst - y).^2, 2)));
    if testPlot
        fprintf(1, 'iter = %f, MVE = %f, alpha = %f\n', iter, MVE, alpha);
    end
   
    %Update Stall Count
    if model.MVE - MVE < 1E-3
        stallCount = stallCount + 1;
        if stallCount > maxStallCount
            break;
        end
    else
        stallCount = 0;
    end
    
    %Update Model if improved
    if MVE < model.MVE
        model.MVE = MVE;
        model.W1 = W1;
        model.W2 = W2;
    else
        W1 = 0.5*W1 + 0.5*model.W1;
        W2 = 0.5*W2 + 0.5*model.W2;
    end      
end

%% Compute Fit Statistics
[~, X] = PredictANN(model.W1, model.W2, x);
[lambda, model.MVE] = BoundLineSearch(@(L) PCVE(X, y, L), 0, 0.5*numSamples, testPlot);
I = eye(model.numHLN + 1);
model.W2 = (pinv(X'*X + lambda*I)*X'*y)';
% yEst = PredictANN(model.W1, model.W2, x);

%Test Plot if Testing Algorithm
if (testPlot)
    if model.InputDim == 1
        figure();
        title('Test Data');
        plot(x, y, '.k');
        hold on;
        x = (min(x):0.001:max(x))';
        yEst = PredictANN(model.W1, model.W2, x);
        plot(x, yEst, '.g');
        
    elseif model.InputDim == 2
        figure();
        PlotModel(model);
        hold on;
        plot3(x(:,1), x(:,2), y, '.k');
        zlim([-3, 3]);
    end
end

end

function err = Error(W1, W2, G1, G2, x, y, alpha)

W1 = W1 - alpha*G1;
W2 = W2 - alpha*G2;
yEst = PredictANN(W1, W2, x);
err = sqrt(mean(sum((yEst - y).^2, 2)));

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
    err = 4*(yEst - y(i,:))^3;
    G2 = G2 + err*g';
    dedg = repmat(err*W2, m, 1)';
    grad = dedg.*((1 - g.^2)*x(i, :));
    G1 = G1 + grad(1:end-1,:);
end

%  norm = 1E-2*sign(W1);
%  norm(:,end) = 0;
G1 = G1/n;% - norm;
G2 = G2/n;    
end


function err = PCVE(X, y, lambda)

%Select a spread of p% of all rows for cross validation
n = size(X,1);
p = min(40/n, 1.0);
s = ceil(1/p);
cvSet = (s:s:n)';
m = length(cvSet);
row = (1:n)';
errVect = zeros(m,1);

%Ensure Matrix is Semi-Positive Definite
lambda = max(lambda, 1E-3);
I = eye(size(X,2));
XtX = X'*X;
[~, p] = chol(XtX + lambda*I);
while p ~= 0
    lambda = lambda*2;
    [~, p] = chol(XtX + lambda*I);
end

%Loop and fit model m times computing error each time
for i = 1:m
   indx = row ~= cvSet(i);
   Xfit = X(indx,:);
   yfit = y(indx);
   temp = inv(chol(Xfit'*Xfit + lambda*I));
   C = (temp*temp')*(Xfit'*yfit);
   yEst = X(~indx,:)*C;
   errVect(i) = (yEst - y(~indx))^2;
end
err = sqrt(mean(errVect));

end



