


function model = RBF(x, y, varargin)
%% Input Processing
if nargin < 2
    error('x and y are required inputs');
elseif size(x,1) ~= size(y,1)
    error('Fitting Inputs and Outputs must have the same size');
end

%% Initialize Values
model = SurrogateModel();
model.Type = 'RBF';
model.InputDim = size(x,2);
model.OutputDim = size(y,2);

numSamples = size(x,1);
crossvalidate = false;
validationSet = false;
xMVE = ones(1, model.InputDim);
yMVE = ones(1, model.InputDim);
testPlot = false;

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
       validationSet = true;
       xMVE = varargin{i+1}{1};
       xMVE = NormalizeInput(xMVE);
       yMVE = varargin{i+1}{2};
       i = i + 3;
       if size(xMVE,2) ~= model.InputDim
           error('Validation set must have same dimension as construction set');
       elseif any(size(xMVE) ~= size(yMVE))
           error('Validation size inconsistent');
       end
   elseif strcmp(varargin{i}, 'semivariance')
       crossvalidate = true;
       i = i + 1;
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
   elseif strcmp(varargin{i}, 'plot')
       testPlot = true;
       i = i + 1;
   else
       error('Unknown command %s', varargin{i});
   end    
end

%% Model Construction Loop
%Generate Gram Matrix
X = zeros(numSamples, numSamples);
dist = zeros(numSamples*(numSamples-1)/2, 1);
vari = zeros(size(dist));
iCount = 0;
for i = 1:numSamples
    for j = i+1:numSamples
        r = max(norm(x(i,:) - x(j,:)), 1E-5);
        iCount = iCount + 1;
        dist(iCount) = r;
        vari(iCount) = 0.5*(y(i) - y(j))^2;
        X(i, j) = r^2*log(r);
        X(j, i) = X(i, j);
    end
end
        
%Determine Regularization Coefficient
if crossvalidate
    %Optimize Regularization Coefficient via Cross Validation
    [model.lambda, model.MVE] = BoundLineSearch(@(L) PCVE(X, y, L), 0, 1, true); 
elseif validationSet
    %Optimize Regularization Coefficient
    [model.lambda, model.MVE] = BoundLineSearch(@(L) VE(X, y, XMVE, yMVE, L), 0, 1);
else       
    %Estimate Lambda from SemiVariance
    bins = ceil(sqrt(numSamples));
    binWidth = max(dist)/bins;
    indx = dist < binWidth;
    while sum(indx) < 3
        binWidth = binWidth*2;
        indx = dist < binWidth;
    end
    R = [ones(sum(indx),1), dist(indx)];
    temp = chol(R'*R)\eye(2);
    fit = (temp*temp')*(R'*vari(indx));
    model.lambda = max(fit(1), 0);
    model.MVE = PCVE(X, y, model.lambda);
end

%Compute Model Parameters and Fit Statistics
I = eye(size(X,2));
temp = chol(X'*X + model.lambda*I)\I;
C = (temp*temp')*(X'*y);
yEst = X*C;
yBar = mean(y);
SStot = sum((y - yBar).^2);
err = yEst - y;
SSres = sum(err.*err);
model.MRSE = sqrt(SSres/length(err));
model.R2 = 1 - SSres/SStot;
model.C = C;
model.nodeX = x;
model.numNodes = numSamples;

%Test Plot if Testing Algorithm
if (testPlot)
    if model.InputDim == 1
        figure();
        plot(x*model.InputSigma + model.InputMu, y*model.OutputSigma + model.OutputMu, '.k');
        hold on;
        plot(x*model.InputSigma + model.InputMu, yEst*model.OutputSigma + model.OutputMu, '.g');            
    elseif model.InputDim == 2
        figure();
        PlotModel(model);
        hold on;
        plot3(x(:,1), x(:,2), y, '.k');
        zlim([-3, 3]);
    end
end
        
end

%Function to compute Partial Cross Validation Error for RBF Regression
function err = PCVE(X,y,lambda)

%Select a spread of p% of all rows for cross validation
n = size(X,1);
p = min(50/n, 1.0);
s = ceil(1/p);
cvSet = (s:s:n)';
m = length(cvSet);
row = (1:n)';
errVect = zeros(m,1);
I = eye(size(X, 2) - 1);

%Update Lambda to ensure positive definiteness
L = max(lambda, 1E-3)*eye(size(X,2));
XtX = X'*X;
[cholX, p] = chol(XtX + L);
while p ~= 0
    L = 2*L;
    [cholX, p] = chol(XtX + L);
end  

%Loop and fit model m times computing error each time
for i = 1:m
   indx = row ~= cvSet(i);
   Xfit = cholX(indx,indx);
   yfit = y(indx);
   temp = Xfit\I;
   C = (temp*temp')*(Xfit'*yfit);
   yEst = X(~indx,indx)*C;
   errVect(i) = (yEst - y(~indx))^2;
end
err = sqrt(mean(errVect));
end

%Compute Validation Error of a Model
function err = VE(X, y, Xval, yval, lambda)

%Update Lambda to ensure positive definiteness
I = eye(size(X,2));
L = max(lambda, 1E-3)*I;
XtX = X'*X;
[cholX, p] = chol(XtX + L);
while p ~= 0
    L = 2*L;
    [cholX, p] = chol(XtX + L);
end  

temp = cholX\I;
C = (temp*temp')*(X*y);
yEst = Xval*C;
err = norm(yEst - yval);

end

