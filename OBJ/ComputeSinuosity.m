function value = ComputeSinuosity(func, varMin, varMax, samplePoints, replications)
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
    samplePoints = 100;
end
derivEp = 1E-4*(varMax - varMin);
dim = length(varMin);
sample = zeros(1,dim);
sampleLoc = rand(samplePoints,dim);
for i = 1:dim
    range = varMax(i) - varMin(i);
    sampleLoc(:,i) = sampleLoc(:,i)*range + varMin(i);
    sample(i) = varMin(i);
end



%% Compute Second Derivative at Sample Locations
d2fdx2 = zeros(samplePoints, dim);
fVals = zeros(3,1);
sumVals = 0;
for i = 1:size(d2fdx2,1)
   for j = 1:dim
       sample = sampleLoc(i,:);
       step = derivEp(j);
       fVals(1) = func(sample);
       sumVals = sumVals + abs(fVals(1));
       sample(j) = sampleLoc(i,j) + step;
       fVals(2) = func(sample);
       sample(j) = sampleLoc(i,j) - step;
       fVals(3) = func(sample);
       d2fdx2(i,j) = (fVals(2) - 2*fVals(1) + fVals(3))/step^2;
   end    
end


%% Compute Sinuosity
meanfVal = sumVals/samplePoints;
mean2ndDeriv = mean(d2fdx2,1);
std2ndDeriv = std(d2fdx2, 0, 1);
value = norm(std2ndDeriv ./ mean2ndDeriv);
% value = norm(var2ndDeriv./meanfVal);

% plot3(sampleLoc(:,1), sampleLoc(:,2), d2fdx2(:,1), '.');


end