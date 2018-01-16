
function yEst = EvalModel(model, x)

%% Error Checking
if ~isstruct(model)
    error('model must be a structure containing surrogate model parameters');
elseif ~isfield(model, 'InputDim') || ~isfield(model, 'OutputDim') || ~isfield(model, 'Type') || ~isfield(model, 'InputMu') || ~isfield(model, 'InputSigma') || ~isfield(model, 'OutputMu') || ~isfield(model, 'OutputSigma')
    error('model should be a product of the SurrogateModel Function');
elseif size(x,2) ~= model.InputDim
    error('x must be an n-by-m matrix where n is the number of points and m is the model input dimension');
end

% Normalize Input
x = NormalizeData(x, model.InputMu, model.InputSigma);
numSamples = size(x,1);

if any(abs(x(:)) > 1.02)
    warning('Extrapolating Surrogate Model.  Results May Be Inaccurate.');
end

%% Model Type Switch
switch model.Type
    case 'PolyReg'
        numCoef = size(model.P,1);
        X = ones(numSamples, numCoef);
        for i = 2:numCoef
            for j = 1:model.InputDim
                X(:,i) = X(:,i).*(x(:,j).^model.P(i,j));
            end
        end
        yEst = X*model.C;
       
    case 'SplineReg'
        numPolyCoef = size(model.P,1);
        numKnotCoef = size(model.K,1)*model.InputDim;
        polyOrder = max(model.P(:));
        X = ones(numSamples, numPolyCoef + numKnotCoef);
        for i = 2:numPolyCoef
            for j = 1:model.InputDim
                X(:,i) = X(:,i).*(x(:,j).^model.P(i,j));
            end
        end
        for i = 1:numKnotCoef
            knot = ceil(i/model.InputDim);
            dim = mod(i-1, model.InputDim) + 1;
            X(:,numPolyCoef+i) = max(x(:,dim) - model.K(knot, dim), 0).^polyOrder;
        end
        yEst = X*model.C;
        
    case 'RBF'
        X = zeros(numSamples, model.numNodes);
        for i = 1:model.numNodes
            nodeVal = repmat(model.nodeX(i,:), numSamples, 1);
            r = sqrt(sum((nodeVal - x).^2, 2));
            r(r < 1E-5) = 1E-5;
            X(:, i) = r.*r.*log(r);
        end
        yEst = X*model.C;
        
    case 'Kriging'
        X = zeros(model.numNodes, numSamples);
        C = -3/(model.Range*model.Range);
        for i = 1:model.numNodes
            nodeVal = repmat(model.nodeX(i,:), numSamples, 1);
            lag = sqrt(sum((nodeVal - x).^2, 2));
            X(i, :) = (model.Sill-model.Nugget)*exp(C*lag.^2) + model.Nugget;
        end
        yEst = (model.nodeY'*model.PsiInv*X)';
        
    case 'ANN'
        yEst = PredictANN(model.W1, model.W2, x);
        
    otherwise
        error('Unsupported Model Type');
end

%Un Studentize Prediction
yEst = yEst*model.OutputSigma + model.OutputMu;

end