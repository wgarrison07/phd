
function model = Kriging(x, y, varargin)
%% Input Processing
if nargin < 2
    error('x and y are required inputs');
elseif size(x,1) ~= size(y,1)
    error('Fitting Inputs and Outputs must have the same size');
end

%% Initialize Values
model = SurrogateModel();
model.Type = 'Kriging';
model.InputDim = size(x,2);
model.OutputDim = size(y,2);

numSamples = size(x,1);
testPlot = false;

%Normalize and Sort Inputs, Normalize Outputs
[x, model.InputMu, model.InputSigma] = NormalizeData(x);
[y, model.OutputMu, model.OutputSigma] = NormalizeData(y);
sumX = sum(x,2);
[~, indx] = sort(sumX);
x = x(indx,:);
y = y(indx);

%% Compute Semi-Variogram
%Compute the lag and semi-variance
gamma = zeros(0.5*numSamples*(numSamples-1),1);
r = zeros(size(gamma));
id = 1;
for i = 1:numSamples
    for j = i+1:numSamples
        gamma(id) = (y(i) - y(j))^2;
        r(id) = norm(x(i, :) - x(j,:));
        id = id + 1;
    end
end

%% Build Covariance Model
bins = ceil(0.5*numSamples);
maxR = 0.99*max(r);
binWidth = maxR/bins;
indx = r < binWidth;
while sum(indx) < 10
    binWidth = binWidth*2;
    indx = r < binWidth;
end
N = median(gamma(indx));
R = BoundLineSearch(@(R) FitSemiVariogram(N, R, r, gamma, 'gaussian'), 0.1*maxR, maxR);
S = mean(gamma(r >= R));
N = min(N, S);

%% Compute Kriging Fit
lag = zeros(numSamples, 1);
Psi = zeros(numSamples, numSamples);
gaussC = -3/(R*R);
for i = 1:numSamples
    for j = 1:numSamples
        lag(j) = norm(x(i,:) - x(j,:));
    end
    Psi(i,:) = (S-N)*exp(gaussC*lag.^2) + N;
end
temp = chol(Psi + N*eye(numSamples))\eye(numSamples);
model.PsiInv = temp*temp';
model.Nugget = N;
model.Range = R;
model.Sill = S;

%% Compute Fit Statistics
yEst = zeros(size(y));
for i = 1:length(yEst)
    yEst(i) = y'*model.PsiInv*Psi(:,i);
end
yBar = mean(y);
SStot = sum((y - yBar).^2);
err = yEst - y;
SSres = sum(err.*err);
model.MRSE = sqrt(SSres/length(err));
model.R2 = 1 - SSres/SStot;
model.numNodes = numSamples;
model.nodeX    = x;
model.nodeY    = y;

%% Plotting
if (testPlot)
    if model.InputDim == 1
        figure();
        plot(x, y, '.k');
        hold on;
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

function [PsiInv, Psi] = ComputeSemiVarianceMatrix(Nugget, Range, Sill, Type, x)

numSamples = size(x,1);
Psi = zeros(numSamples, numSamples);
N = Nugget;
S = Sill;
R = Range;
for i = 1:numSamples
    for j = i:numSamples
        lag = norm(x(i,:) - x(j,:));
        Psi(i, j) = S - VarianceModel(N, S, R, lag, Type);
        Psi(j, i) = Psi(i, j);
    end
end
 temp = chol(Psi + N*eye(numSamples))\eye(numSamples);
 PsiInv = temp*temp';
 
end

function [err] = FitSemiVariogram(Nugget, Range, r, gamma, Type)

%Compute Sill from data
indx = r >= Range;
if sum(indx) == 0
    Sill = gamma(end);
else
    Sill = mean(gamma(indx));
end

%Compute Fit Error
err = (VarianceModel(Nugget, Sill, Range, r, Type) - gamma).^2;
err = sqrt(mean(err));

end

