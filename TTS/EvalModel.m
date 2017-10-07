
function yEst = EvalModel(model, x)

%% Error Checking
if ~isstruct(model)
    error('model must be a structure containing surrogate model parameters');
elseif ~isfield(model, 'InputDim') || ~isfield(model, 'OutputDim') || ~isfield(model, 'Type') || ~isfield(model, 'Mu') || ~isfield(model, 'Sigma')
    error('model should be a product of the SurrogateModel Function');
elseif size(x,2) ~= model.InputDim
    error('x must be an n-by-m matrix where n is the number of points and m is the model input dimension');
end

% Normalize Input
x = NormalizeInput(x, model.Mu, model.Sigma);
numSamples = size(x,1);

%% Model Type Switch
switch model.Type
    case 'PolyReg'
        %% Polynomial Regression Model
        numCoef = size(model.P,1);
        X = ones(numSamples, numCoef);
        for i = 2:numCoef
            for j = 1:model.InputDim
                X(:,i) = X(:,i).*(x(:,j).^model.P(i,j));
            end
        end
        yEst = X*model.C;
       
    case 'SplineReg'
        %% Spline Regression Model
        numPolyCoef = size(model.P,1);
        numKnotCoef = size(model.K,1);
        polyOrder = max(model.P(:));
        X = ones(numSamples, numPolyCoef + numKnotCoef);
        for i = 2:numPolyCoef
            for j = 1:model.InputDim
                X(:,i) = X(:,i).*(x(:,j).^model.P(i,j));
            end
        end
        for i = 1:numKnotCoef
            X(:,i+numPolyCoef) = sum(max(x - model.K(i), 0),2).^polyOrder;
        end
        yEst = X*model.C;
        
    otherwise
    %% Unsupported Model Type
        error('Unsupported Model Type');
end


    

end