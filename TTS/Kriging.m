
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
crossValidate = true;
xMVE = ones(1, model.InputDim);
yMVE = ones(1, model.InputDim);
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
    sum(indx)
    gammaBar(i) = mean(gamma(indx));
end
binR = binR(1:end-1);
gammaBar = gammaBar(1:end-1);
fprintf(1, 'Bin 1 SNR = %f\n', 1/gammaBar(1));

%Compute the Range (Sill is Regressed)
R = binR(end);
for i = 1:length(binR)-1
    if gammaBar(i+1) < gammaBar(i)
        R = mean(binR(i:end));
        break;
    end
end
plot(binR, gammaBar, '.k');


%Determine whether to use a Spherical, Exponential, or Gaussian empirical
%semi-variogram model.
%Spherical Model
nLag = binR/R;
sphericalV = 3/2*min(nLag,1) - 0.5*min(nLag,1).^3;
X = [sphericalV, 1-sphericalV];
sphericalC = (X'*X)\(X'*gammaBar);
sphericalC(2) = max(sphericalC(2), 0);
sphericalGamma = X*sphericalC;
sphericalErr = norm(sphericalGamma - gammaBar);
hold on;
plot(binR, sphericalGamma, '-g');

%Exponential Model
exponentialV = 1 - exp(-3*nLag);
X = [exponentialV, 1-exponentialV];
exponentialC = (X'*X)\(X'*gammaBar);
exponentialC(2) = max(exponentialC(2), 0);
exponentialGamma = X*exponentialC;
exponentialErr = norm(exponentialGamma - gammaBar);
plot(binR, exponentialGamma, '-b');

%Gaussian Model
gaussianV = 1 - exp(-3*nLag.^2);
X = [gaussianV, 1-gaussianV];
gaussianC = (X'*X)\(X'*gammaBar);
gaussianC(2) = max(gaussianC(2), 0);
gaussianGamma = X*gaussianC;
gaussianErr = norm(gaussianGamma - gammaBar);
plot(binR, gaussianGamma, '-r');

%Select the Variogram Type
model.Range = R;
if sphericalErr < exponentialErr && sphericalErr < gaussianErr
    model.vModel = 'spherical';
    model.Sill   = sphericalC(1);
    model.Nugget = sphericalC(2);
elseif exponentialErr < sphericalErr && exponentialErr < gaussianErr
    model.vModel = 'exponential';
    model.Sill   = exponentialC(1);
    model.Nugget = exponentialC(2);
else
    model.vModel = 'gaussian';
    model.Sill   = gaussianC(1);
    model.Nugget = gaussianC(2);
end

%% Compute Model Parameters
model.Psi = zeros(numSamples, numSamples);
N = model.Nugget;
fprintf(1, 'Estimated Data SNR = %f\n', 1/N);
S = model.Sill;
R = model.Range;
for i = 1:numSamples
    for j = i:numSamples
        lag = min(norm(x(i,:) - x(j,:))/R, 1);
        if     strcmp(model.vModel, 'spherical')
            model.Psi(i, j) = (S-N)*(3/2*lag - 0.5*lag^3);
            model.Psi(j, i) = model.Psi(i, j);
        elseif strcmp(model.vModel, 'exponential')
            model.Psi(i, j) = (S-N)*(1 - exp(-3*lag));
            model.Psi(j, i) = model.Psi(i, j);
        elseif strcmp(model.vModel, 'gaussian')
            model.Psi(i, j) = (S-N)*(1 - exp(-3*lag^2));
            model.Psi(j, i) = model.Psi(i, j);
        else
            error('Unsupported Model Type');
        end
    end
end
model.PsiInv = inv(model.Psi - N*eye(numSamples));

%% Compute Fit Statistics
yEst = zeros(size(y));
for i = 1:length(yEst)
    yEst(i) = y'*model.PsiInv*(model.Psi(:,i));
end
figure();
plot(x, y, '.k');
hold on;
plot(x, yEst, '.g');


end

function err = FitSemiVariogram(Nugget, Type)


end