
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
            for j = 1:numSamples
                r = max(norm(model.nodeX(i,:) - x(j,:)), 1E-5);
                X(j, i) = r^2*log(r);
            end
        end
        yEst = X*model.C;
        
    case 'Kriging'
        X = zeros(model.numNodes, numSamples);
        for i = 1:model.numNodes
            for j = 1:numSamples
                lag = min(norm(model.nodeX(i,:) - x(j,:))/model.Range, 1);
                s = VarianceModel(model.Nugget, model.Sill, lag, model.vModel);
                X(i, j) = s;
            end
        end
        yEst = model.nodeY'*model.PsiInv*X;
        
    case 'ANN'
        yEst = PredictANN(model.W1, model.W2, x);
        
    otherwise
        error('Unsupported Model Type');
end

%Un Studentize Prediction
yEst = yEst*model.OutputSigma + model.OutputMu;

end