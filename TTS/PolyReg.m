

function model = PolyReg(x, y, varargin)
%% Input Processing
if nargin < 2
    error('x and y are required inputs');
elseif size(x,1) ~= size(y,1)
    error('Fitting Inputs and Outputs must have the same size');
end

%% Initialize Values
model = SurrogateModel();
model.Type = 'PolyReg';
model.InputDim = size(x,2);
model.OutputDim = size(y,2);
model.MVE   = 1E16;
model.C     = 0;
model.P     = zeros(1, model.InputDim);

numSamples = size(x,1);
crossValidate = true;
xMVE = ones(1, model.InputDim);
yMVE = ones(1, model.InputDim);
startOrder = 0;
stopOrder = 0;
while ComputePowerMatrixSize(stopOrder + 1, model.InputDim) < numSamples
    stopOrder = stopOrder + 1;
end
testPlot = false;
abortTol = 1E-3;
maxStallCount = 3;
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
   elseif strcmp(varargin{i}, 'plot')
       testPlot = true;
       i = i + 1;
   else
       error('Unknown command %s', varargin{i});
   end    
end

%% Model Construction Loop
for iOrder  = startOrder:stopOrder
    %Generate Vandermoode Matrix
    P = ComputePowerMatrix(iOrder, model.InputDim);
    numCoef = size(P,1);
    X = ones(numSamples, numCoef);
    XMVE = ones(size(xMVE,1), numCoef);
    for i = 2:numCoef
        for j = 1:model.InputDim
            X(:,i) = X(:,i).*(x(:,j).^P(i,j));
            XMVE(:,i) = XMVE(:,i).*(xMVE(:,j).^P(i,j));
        end
    end

    try %Fit Model
        if crossValidate
            %Optimize Regularization Coefficient via Cross Validation
            [MVE, C] = PCVE(X, y); 
        else       
            %Optimize Regularization Coefficient
            [MVE, C] = VE(X, y, XMVE, yMVE);         
        end
    catch %If Cholesky fails, break loop
        break;
    end

    %Compute Model Parameters and Fit Statistics
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
            tempModel.C = C;
            tempModel.P = P;
            
            figure();
            PlotModel(tempModel);
            hold on;
            plot3(x(:,1), x(:,2), y, '.k');
            zlim([-1, 1]);
        end
    end
        
    %Update Model Parameters if Error Decreased
    if MVE < model.MVE
        model.Order = iOrder;
        model.C = C;
        model.P = P;
        model.MRSE = MRSE;
        model.R2 = R2;
        model.MVE = MVE;
        stallCount = 0;
    else
        stallCount = stallCount + 1;
    end
        
    if testPlot
        fprintf(1, 'Order %i, MVE = %f\n', iOrder, MVE);
    end
    
    %Break generation loop if if model meets abort criteria
    if model.MVE < abortTol || stallCount > maxStallCount
        break;
    end
end

end

%Function to compute Partial Cross Validation Error for Polynomial
%Regression.
function [err, Cavg] = PCVE(X,y)

%Select a spread of p% of all rows for cross validation
n = size(X,1);
p = min(40/n, 1.0);
s = ceil(1/p);
cvSet = (s:s:n)';
m = length(cvSet);
row = (1:n)';
errVect = zeros(m,1);
I = eye(size(X,2));
Cavg = zeros(size(X,2),1);

%Loop and fit model m times computing error each time
for i = 1:m
   indx = row ~= cvSet(i);
   Xfit = X(indx,:);
   yfit = y(indx);
   temp = chol(Xfit'*Xfit)\I;
   C = (temp*temp')*(Xfit'*yfit);
   Cavg = Cavg + C;
   yEst = X(~indx,:)*C;
   errVect(i) = (yEst - y(~indx))^2;
end
err = sqrt(mean(errVect));
Cavg = Cavg/m;
end

%Compute Validation Error of a Model
function [err, C] = VE(X, y, Xval, yval)

temp = chol(X'*X)\eye(size(X,2));
C = (temp*temp')*(X*y);
yEst = Xval*C;
err = norm(yEst - yval);

end



