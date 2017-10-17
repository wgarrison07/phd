
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
        r(id) = norm(x(i, :) - x(j,:));
        id = id + 1;
    end
end

%Bin Data
range = max(r);
bins = ceil(sqrt(numSamples));
binR = (0:range/bins:range)';
gammaBar = zeros(length(binR),1);
for i = 1:length(gammaBar)-1
    indx = r >= binR(i) & r < binR(i+1);
    gammaBar(i) = mean(gamma(indx));
end
binR = binR(1:end-1);
gammaBar = gammaBar(1:end-1);

%Compute the Range and Sill
R = binR(end);
S = gammaBar(end);
for i = 1:length(binR)-1
    if gammaBar(i+1) < gammaBar(i)
        R = binR(i);
        S = mean(gammaBar(i:end));
        break;
    end
end

%Determine whether to use a Spherical, Exponential, or Gaussian empirical
%semi-variogram model.
%Spherical Model
sphericalErr = 1E9;
% [sphericalN, sphericalErr] = fminbnd(@(N) FitSemiVariogram(N, R, S, x, y, 'spherical'), 0, S);

%Exponential Model
exponentialErr = 1E9;
% [exponentialN, exponentialErr] = fminbnd(@(N) FitSemiVariogram(N, R, S, x, y, 'exponential'), 0, S);

%Gaussian Model
[gaussianN, gaussianErr] = fminbnd(@(N) FitSemiVariogram(N, R, S, x, y, 'gaussian'), 0, S, optimset('Display', 'iter'));

%Select the Variogram Type
model.Range = R;
model.Sill  = S;
if sphericalErr < exponentialErr && sphericalErr < gaussianErr
    model.vModel = 'spherical';
    model.Nugget = sphericalN;
    model.MVE    = sphericalErr;
elseif exponentialErr < sphericalErr && exponentialErr < gaussianErr
    model.vModel = 'exponential';
    model.Nugget = exponentialN;
    model.MVE    = exponentialErr;
else
    model.vModel = 'gaussian';
    model.Nugget = gaussianN;
    model.MVE    = gaussianErr;
end

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

function [PsiInv, Psi] = ComputeSemiVarianceMatrix(Nugget, Range, Sill, Type, x)

numSamples = size(x,1);
Psi = zeros(numSamples, numSamples);
N = Nugget;
S = Sill;
R = Range;
for i = 1:numSamples
    for j = i:numSamples
        lag = min(norm(x(i,:) - x(j,:))/R, 1);
        if  strcmp(Type, 'spherical')
            Psi(i, j) = (S-N)*(3/2*lag - 0.5*lag^3);
            Psi(j, i) = Psi(i, j);
        elseif strcmp(Type, 'exponential')
            Psi(i, j) = (S-N)*(1 - exp(-3*lag));
            Psi(j, i) = Psi(i, j);
        elseif strcmp(Type, 'gaussian')
            Psi(i, j) = (S-N)*(1 - exp(-3*lag^2));
            Psi(j, i) = Psi(i, j);
        else
            error('Unsupported Model Type');
        end
    end
end
PsiInv = inv(Psi - N*eye(numSamples));

end

function [err] = FitSemiVariogram(Nugget, Range, Sill, x, y, Type)

%Compute covariance matrix and inverse
[PsiInv, Psi] = ComputeSemiVarianceMatrix(Nugget, Range, Sill, Type, x);

%Compute Crossvalidation Error
numSamples = length(y);
errVect = zeros(numSamples, 1);
iVect = 1:numSamples;
for i = iVect
    indx = iVect ~= i;
    Pinv = PsiInv(indx, indx);
    yEst = y(indx)'*Pinv*Psi(indx, i);
    errVect(i) = (yEst - y(i))^2;
end
err = sqrt(mean(errVect));

end