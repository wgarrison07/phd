

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

%If less than 3rd order data points, just build a polynomial regression
polyOrder = 3;
numSamples = size(x,1);
crossValidate = true;
xMVE = ones(1, model.InputDim);
yMVE = ones(1, model.InputDim);
startOrder = 1;
stopOrder = (numSamples - 1) - polyOrder^model.InputDim;
testPlot = true;
abortTol = 1E-3;
maxStallCount = 20;
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
       if startOrder >= stopOrder
           error('startOrder larger than degrees of freedom');
       end
       i = i + 2;
   elseif strcmp(varargin{i}, 'order') 
       startOrder = varargin{i+1};
       if (startOrder > stopOrder)
           error('order larger than degrees of freedom');
       else
           stopOrder = startOrder;
       end
       i = i + 2;
   elseif strcmp(varargin{i}, 'polyOrder')
       polyOrder = varargin{i+1};
       stopOrder = (numSamples - 1) - polyOrder^model.InputDim;
       i = i + 2;
   else
       error('Unknown command %s', varargin{i});
   end    
end

%% If Insuffcient Number of Samples, Use Polynomial Regression
if (numSamples - 2 < polyOrder^model.InputDim)
    model = PolyReg(x, y);
    model.Type = 'SplineReg';
    return;
end

%% Model Construction Loop
model.MVE   = Inf;
for iOrder  = startOrder:stopOrder
    % Generate Knot Positions
    step = 1/(iOrder+1);
    knotIndx = ceil(numSamples*(step:step:(1-step)));
    K = x(knotIndx,:);
    
    % Compute Vandermoode Matrix
    P = ComputePowerMatrix(polyOrder, model.InputDim);
    numPolyCoef = size(P,1);
    numKnotCoef = iOrder;
    X = ones(numSamples, numPolyCoef + numKnotCoef);
    for i = 2:numPolyCoef
        for j = 1:model.InputDim
            X(:,i) = X(:,i).*(x(:,j).^P(i,j));
        end
    end
    for i = 1:numKnotCoef
        X(:,i+numPolyCoef) = sum(max(x - K(i), 0),2).^polyOrder;
    end
    
    if crossValidate
        %Optimize Regularization Coefficient via Cross Validation
        [lambda, MVE] = fminbnd(@(L) PCVE(X, y, L), 0, numSamples); 
    else       
        %Optimize Regularization Coefficient
        [lambda, MVE] = fminbnd(@(L) VE(X, y, XMVE, yMVE, L), 0, numSamples);         
    end

    %Compute Model Parameters and Fit Statistics
    L = lambda*eye(size(X,2));
    C = (X'*X + L)\(X'*y);
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
            figure();
            plot3(x(:,1), x(:,2), y, '.k');
            hold on;
            plot3(x(:,1), x(:,2), yEst, '.g');
        end
    end
        
    %Update Model Parameters if Error Decreased
    fprintf(1, 'Order %i, MVE = %f\n', iOrder, MVE);
    if MVE < model.MVE
        model.Order = iOrder;
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

    %Break generation loop if if model meets abort criteria
    if model.MVE < abortTol || stallCount > maxStallCount
        break;
    end
end
    
end

%Function to compute Partial Cross Validation Error for Polynomial
%Regression.
function err = PCVE(X,y,lambda)

%Select a spread of p% of all rows for cross validation
p = 1;
n = size(X,1);
s = ceil(1/p);
cvSet = (s:s:n)';
m = length(cvSet);
row = (1:n)';
errVect = zeros(m,1);
L = lambda*eye(size(X,2));

%Loop and fit model m times computing error each time
for i = 1:m
   indx = row ~= cvSet(i);
   Xfit = X(indx,:);
   yfit = y(indx);
   C = (Xfit'*Xfit + L)\(Xfit'*yfit);
   yEst = X(~indx,:)*C;
   errVect(i) = yEst - y(~indx);
end
err = norm(errVect/m);

end

%Compute Validation Error of a Model
function err = VE(X, y, Xval, yval, lambda)

L = lambda*eye(size(X,2));
C = (X'*X + L)\(X*y);
yEst = Xval*C;
err = norm(yEst - yval);

end
