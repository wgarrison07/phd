
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
testPlot = true;

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
        gamma(id) = 0.5*(y(i) - y(j))^2;
        r(id) = 0.5*norm(x(i, :) - x(j,:));
        id = id + 1;
    end
end

%Bin SemiVariance and Compute Nugget
bins = ceil(0.5*numSamples);
binR = zeros(bins,1);
gammaBar = zeros(bins, 1);
binWidth = max(r)/bins;
for i = 1:bins
    indx = r < i*binWidth & r >= (i-1)*binWidth;
    binR(i) = mean(r(indx));
    gammaBar(i) = mean(gamma(indx));
end
N = gammaBar(1); %Temporary Nugget

%% Determine whether to use a Linear, Spherical, Exponential, or Gaussian Model
[linearR, linearErr] = BoundLineSearch(@(R) FitSemiVariogram(N, R, binR, gammaBar, 'linear'), 0, binR(end));
[sphericalR, sphericalErr] = BoundLineSearch(@(R) FitSemiVariogram(N, R, binR, gammaBar, 'spherical'), 0, binR(end));
[exponentialR, exponentialErr] = BoundLineSearch(@(R) FitSemiVariogram(N, R, binR, gammaBar, 'exponential'), 0, binR(end));
[gaussianR, gaussianErr] = BoundLineSearch(@(R) FitSemiVariogram(N, R, binR, gammaBar, 'gaussian'), 0, binR(end));

%Select the Variogram Type
if linearErr < sphericalErr && linearErr < exponentialErr && linearErr < gaussianErr
    model.vModel = 'linear';
    model.Range  = linearR;
elseif  sphericalErr < linearErr && sphericalErr < exponentialErr && sphericalErr < gaussianErr
    model.vModel = 'spherical';
    model.Range  = sphericalR;
elseif exponentialErr < linearErr && exponentialErr < sphericalErr && exponentialErr < gaussianErr
    model.vModel = 'exponential';
    model.Range  = exponentialR;
else
    model.vModel = 'gaussian';
    model.Range = gaussianR;
end
if any(r >= model.Range)
    model.Sill = mean(gamma(r >= model.Range));
else
    model.Sill = gammaBar(end);
end

%% Compute the Regression Nugget via CrossValidation
[model.Nugget, model.MVE] = BoundLineSearch(@(N) PCVE(N, model.Range, model.Sill, model.vModel, x, y), N, model.Sill, true);
% model.Nugget = N;

%% Compute Model Parameters & Fit Statistics
[model.PsiInv, Psi] = ComputeSemiVarianceMatrix(model.Nugget, model.Range, model.Sill, model.vModel, x);
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
        plot3(x(:,1), x(:,2), y, '.k');
        hold on;
        plot3(x(:,1), x(:,2), yEst, '.g');
    end
end


end

function err = PCVE(Nugget, Range, Sill, Type, x, y)
%Select a spread of p% of all rows for cross validation
p = 1;
n = size(x,1);
s = ceil(1/p);
cvSet = (s:s:n)';
m = length(cvSet);
row = (1:n)';
errVect = zeros(m,1);
[~, Psi] = ComputeSemiVarianceMatrix(Nugget, Range, Sill, Type, x);
N = Nugget*eye(size(Psi)-1);
yt = y';

%Loop and fit model m times computing error each time
for i = 1:m
   indx = row ~= cvSet(i);
   subPsi = Psi(indx, indx);
   PsiInv = inv(subPsi - N);
   yEst = yt(indx)*PsiInv*Psi(indx, cvSet(i));
   errVect(i) = (yEst - y(~indx))^2;
end
err = sqrt(mean(errVect));


end

function [PsiInv, Psi] = ComputeSemiVarianceMatrix(Nugget, Range, Sill, Type, x)

numSamples = size(x,1);
Psi = zeros(numSamples, numSamples);
N = min(Nugget, Sill);
S = Sill;
R = Range;
for i = 1:numSamples
    for j = i:numSamples
        lag = min(norm(x(i,:) - x(j,:))/R, 1);
        Psi(i, j) = VarianceModel(N, S, lag, Type);
        Psi(j, i) = Psi(i, j);
    end
end
PsiInv = inv(Psi - N*eye(numSamples));

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
lag = min(r/Range, 1);
err = (VarianceModel(Nugget, Sill, lag, Type) - gamma).^2;
err = sqrt(mean(err));

end

