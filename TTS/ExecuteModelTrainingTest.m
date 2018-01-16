function [R2, samples] = ExecuteModelTrainingTest(ObjFunction, ModelType, SNR, count, totCount, timer)

rng(2194851);
supportedModels = {'PolyReg', 'SplineReg', 'RBF', 'Kriging', 'ANN'};
samples = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000];
reps    = 100;

%% Input Processing
if ~isstruct(ObjFunction)
    if ischar(ObjFunction)
        objs = ObjectiveFunctionStructure();
        for i = 1:length(objs)
            if strcmp(ObjFunction, objs{i}.Name)
                ObjFunction = objs{i};
                break;
            end
        end
        if ischar(ObjFunction)
            error('Unable to find objecive function structure named %s\n', ObjFunction);
        end
    else
        error('ObjFunction must either be a objective function structure or name of one in ObjectiveFunctionStructure()\n');
    end
end

if ~ismember(ModelType, supportedModels)
    error('Only the following ModelType are supported: %s', supportedModels{:});
end

if SNR <= 0
    error('SNR must be > 0');
end

if nargin >= 6
    count = count*reps;
    totCount = totCount*reps;
end

%% Model Training Loop
R2 = zeros(1, length(samples));
maxRows = ceil(samples(end)*1.1);
rows = 1:maxRows;
for i = 1:reps
    %Create Data Set(randomized for each rep
    [inputSet, outputSet] = ProduceDataSet(maxRows, ObjFunction);
    
    %Create Noise Vector
    signal = std(outputSet);
    noise = signal*signal/SNR;
    noiseVals = sqrt(noise)*randn(maxRows, 1);
    
    %Print Status
    if nargin < 6
        fprintf(1, 'Training Rep Number %i\n', i);
    else
        time = toc(timer);
        remainingMin = (totCount - count - i)/(count + i) * (time/60);
        fprintf(1, 'Training Rep Number %i, %5.2f %% Complete, %f min remaining\n', i, (count+i)/totCount*100, remainingMin); 
    end
    order = 0;
    
    for j = 1:length(samples)
        %Determine Sample Index
        indx = rows <= samples(j);
        
        %Construct Model with Noised Data
        model = ConstructModel(ModelType, inputSet(indx,:), outputSet(indx,:) + noiseVals(indx,:), order);
        if ismember(ModelType, supportedModels(1:2))
            order = model.Order;
        end
        
        %Check Model Accuracy against Unnoised Data
        y = outputSet(~indx);
        yEst = EvalModel(model, inputSet(~indx,:));
        yBar = mean(y);
        SStot = sum((y - yBar).^2);
        err = yEst - y;
        SSres = sum(err.*err);
        r2 = 1 - SSres/SStot;
        R2(j) = R2(j) + r2;
        fprintf(1, 'Samples: %f, Err: %f\n', samples(j), r2);
        
    end
    
end

%Normalize Error by Num Replications
R2 = R2/reps;


end

function [inputSet, outputSet] = ProduceDataSet(numSamples, ObjFunction)
    %Determine Sample Range & Construct Inputs
    range = [ObjFunction.MinVals; ObjFunction.MaxVals];
    inputSet = RandomLatinHypercube(range, numSamples);
    
    %Compute outputSet
    outputSet = zeros(numSamples, 1);
    for i = 1:numSamples
        outputSet(i) = ObjFunction.Func(inputSet(i,:));
    end   

end

function model = ConstructModel(type, inputData, outputData, order)
    switch type
        case 'PolyReg'
            model = PolyReg(inputData, outputData, 'startOrder', order);
        case 'SplineReg'
            model = SplineReg(inputData, outputData, 'startOrder', order);
        case 'RBF'
            model = RBF(inputData, outputData);
        case 'Kriging'
            model = Kriging(inputData, outputData);
        case 'ANN'
            model = ANN(inputData, outputData);
        otherwise
            error('Unsupported Model Type %s\n', type);
    end
end
