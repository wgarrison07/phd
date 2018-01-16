function value = ComputeSinuosity(func, varMin, varMax, samplePoints)
%% Check Inputs
if nargin < 3
    error('"func", "varMin", and "varMax" are required inputs');
elseif ~isa(func, 'function_handle')
    error('"func" must be a function handle');
elseif ~all(size(varMin) == size(varMax))
    error('"varMin" and "varMax" must be the same size');
elseif size(varMin,1) > 1 && size(varMin,2) > 1
    error('"varMin" and "varMax" must be 1xn or nx1 vectors that set the variable ranges');
end

%% Initialize Variables
if nargin < 4
    samplePoints = 10000;
end
dim = length(varMin);
sampleLoc = rand(samplePoints,dim);

%% Compute Second Derivative at Sample Locations
d2fdx2 = zeros(samplePoints, dim);
fVals = zeros(3,1);
data = zeros(samplePoints, 1);
stepVar = (varMax - varMin)/1E2;
for i = 1:size(d2fdx2,1)
   for j = 1:dim
       sample = sampleLoc(i,:);
       step = stepVar(j);
       fVals(1) = func(sample);
       sample(j) = sampleLoc(i,j) + step;
       fVals(2) = func(sample);
       sample(j) = sampleLoc(i,j) - step;
       fVals(3) = func(sample);
       d2fdx2(i,j) = (fVals(2) - 2*fVals(1) + fVals(3))/step^2;
       data(i) = fVals(1);
   end    
end


%% Filter Noisy Second Derivative
normDat = NormalizeData(d2fdx2);
indx = all(abs(normDat) < 6, 2);
d2fdx2 = d2fdx2(indx,:);
samplePoints = sum(indx);

%% Determine Optimal Bin Size for Data
% Freedman-Diaconis Rule
range = max(max(d2fdx2) - min(d2fdx2), 1E-2);
innerRange = max(iqr(timeseries(d2fdx2)), 1E-2);
bins = ceil(mean(range ./ (2*innerRange*samplePoints^(-1/3))));
bins = min(bins, ceil(2*sqrt(samplePoints)));
% bins = sqrt(samplePoints);

%% Bin Samples
binIndx = zeros(size(d2fdx2));
for i = 1:dim
    start = min(d2fdx2(:, i));
    width = range(i)/bins;
    for j = 1:bins
        indx = d2fdx2(:, i) > (start + (j-1)*width) & d2fdx2(:, i) <= (start + j*width);
        binIndx(indx, i) = j;
    end
end

%% Compute Joint Entropy
H = 0;
comb = unique(binIndx, 'rows');
for i = 1:size(comb, 1)
    indx = true(samplePoints, 1);
    for j = 1:dim
        indx = indx & (binIndx(:, j) == comb(i, j));
    end
    count = sum(indx);
    p = count/samplePoints;
    H = H + p*log2(p)/log2(samplePoints);
end
value = -1*H;

end