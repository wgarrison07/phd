
function model = SplineReg(x, y, varargin)
%% Input Processing
if nargin < 2
    error('x and y are required inputs');
elseif size(x,1) ~= size(y,1)
    error('Fitting Inputs and Outputs must have the same size');
end

%% Initialize Values
model = SurrogateModel();
model.Type = 'SplineReg';
model.InputDim = size(x,2);
model.OutputDim = size(y,2);

%Truncate Spline Polynomial Order to fit number of samples
numSamples = size(x,1);
polyOrder = 3;
while polyOrder > 0 && (numSamples - 2) < ComputePowerMatrixSize(polyOrder, model.InputDim);
    polyOrder = polyOrder - 1;
end
crossValidate = true;
xMVE = ones(1, model.InputDim);
yMVE = ones(1, model.InputDim);
startOrder = 0;
stopOrder = floor(((numSamples - 1) - ComputePowerMatrixSize(polyOrder, model.InputDim))/model.InputDim);
testPlot = false;
abortTol = 1E-3;
maxStallCount = 5;
stallCount = 0;

%Normalize and Sort Inputs, Normalize Outputs
[x, model.InputMu, model.InputSigma] = NormalizeData(x);
[y, model.OutputMu, model.OutputSigma] = NormalizeData(y);
sumX = sum(x,2);
[~, indx] = sort(sumX);
x = x(indx,:);
y = y(indx);


%% Process Varargin
i = 1;
while i <= length(varargin)
   if strcmp(varargin{i}, 'validationSet')
       crossValidate = false;
       xMVE = varargin{i+1}{1};
       xMVE = NormalizeInput(xMVE);
       yMVE = varargin{i+1}{2};
       i = i + 3;
       if size(xMVE,2) ~= model.InputDim
           error('Validation set must have same dimension as construction set');
       elseif any(size(xMVE) ~= size(yMVE))
           error('Validation size inconsistent');
       end
   elseif strcmp(varargin{i}, 'startOrder')
       startOrder = varargin{i+1};
       if startOrder > stopOrder
           startOrder = 0;
       end
       i = i + 2;
   elseif strcmp(varargin{i}, 'order') 
       startOrder = varargin{i+1};
       if (startOrder > stopOrder)
           startOrder = stopOrder;
       else
           stopOrder = startOrder;
       end
       i = i + 2;
   elseif strcmp(varargin{i}, 'polyOrder')
       polyOrder = varargin{i+1};
       stopOrder = (numSamples - 1) - polyOrder^model.InputDim;
       i = i + 2;
   elseif strcmp(varargin{i}, 'plot')
       testPlot = true;
       i = i + 1;
   else
       error('Unknown command %s', varargin{i});
   end    
end

%% Model Construction Loop
model.MVE   = 1E16;
for iOrder  = startOrder:stopOrder  
    
    % Generate Knot Positions
    step = 1/(iOrder+1);
    knotIndx = ceil(numSamples*(step:step:(1-step)));
    K = x(knotIndx,:);
    
    % Compute Vandermoode Matrix
    P = ComputePowerMatrix(polyOrder, model.InputDim);
    numPolyCoef = size(P,1);
    numKnotCoef = iOrder*model.InputDim;
    X = ones(numSamples, numPolyCoef + numKnotCoef);
    for i = 2:numPolyCoef
        for j = 1:model.InputDim
            X(:,i) = X(:,i).*(x(:,j).^P(i,j));
        end
    end
    for i = 1:numKnotCoef
        knot = ceil(i/model.InputDim);
        dim = mod(i-1, model.InputDim) + 1;
        X(:,numPolyCoef+i) = max(x(:,dim) - K(knot, dim), 0).^polyOrder;
    end
    
    try %Fit Model
        if crossValidate
            [MVE, C] = PCVE(X, y); 
        else       
            [MVE, C] = VE(X, y, XMVE, yMVE);         
        end
    catch
        break;
    end

    %Compute Model Fit Statistics
    yEst = X*C;
    yBar = mean(y);
    SStot = sum((y - yBar).^2);
    err = yEst - y;
    SSres = sum(err.*err);
    MRSE = sqrt(SSres/length(err));
    R2 = 1 - SSres/SStot;
    
    %Test Plot if Testing Algorithm
    if (testPlot)
        if model.InputDim == 1
            figure();
            plot(x, y, '.k');
            hold on;
            plot(x, yEst, '.g');
        elseif model.InputDim == 2
            tempModel = model;
            tempModel.K = K;
            tempModel.C = C;
            tempModel.P = P;
            tempModel.lambda = lambda;            
            
            figure();
            PlotModel(tempModel);
            hold on;
            plot3(x(:,1), x(:,2), y, '.k');
            zlim([-3, 3]);
        end
    end
        
    if model.MVE - MVE > 0.1*model.MVE
        model.Order = iOrder;
        model.K = K;
        model.C = C;
        model.P = P;
        model.lambda = lambda;
        model.MRSE = MRSE;
        model.R2 = R2;
        model.MVE = MVE;
        stallCount = 0;
    else
        stallCount = stallCount + 1;
    end

    %Update Model Parameters if Error Decreased
    if testPlot
        fprintf(1, 'Order %i, MVE = %f, StallCount = %f\n', iOrder, MVE, stallCount);
    end
    
    %Break generation loop if if model meets abort criteria
    if model.MVE < abortTol || stallCount > maxStallCount
        break;
    end
    

end
    
end

%Function to compute Partial Cross Validation Error for Polynomial
%Regression.
function [err, lambda] = PCVE(X,y)

%Select a spread of p% of all rows for cross validation
n = size(X,1);
p = min(50/n, 1.0);
s = ceil(1/p);
cvSet = (s:s:n)';
m = length(cvSet);
row = (1:n)';
errVect = zeros(m,1);

%Ensure Matrix is Positive Definite
I = eye(size(X,2));
L = 1E-3*I;
XtX = X'*X;
[~, p] = chol(XtX + L);
while p ~= 0
    L = 2*L;
    [~, p] = chol(XtX + L);
end
lambda = L(1);

%Loop and fit model m times computing error each time
for i = 1:m
   indx = row ~= cvSet(i);
   Xfit = X(indx,:);
   yfit = y(indx);
   temp = chol(Xfit'*Xfit + L)\I;
   C = (temp*temp')*(Xfit'*yfit);
   yEst = X(~indx,:)*C;
   errVect(i) = (yEst - y(~indx))^2;
end
err = sqrt(mean(errVect));

end

%Compute Validation Error of a Model
function [err, lambda] = VE(X, y, Xval, yval)

%Ensure Matrix is Semi-Positive Definite
lambda = 0;
I = eye(size(X,2));
XtX = X'*X;
[~, p] = chol(XtX + lambda*I);
while p ~= 0
    lambda = lambda + 1;
    [~, p] = chol(XtX + lambda*I);
end

temp = inv(chol(X'*X + lambda*I));
C = (temp*temp')*(X'*y);
yEst = Xval*C;
err = sqrt(mean((yEst - yval).^2));

end
