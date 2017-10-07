%% Objective Evaluation Function
function [valVect, cov] = EvalObjFunc(objFunc, inputs, replications)

    if nargin < 3
        replications = 1;
    elseif replications < 1
        error('replications must be >= 1');
    end
    
    if size(inputs,2) > size(inputs, 1)
        inputs = inputs';
    end
    
    inputDim = size(inputs, 1);
    cStep = 1E-3;
    Xsize = replications*2*inputDim+1;
    X = [ones(Xsize,1), zeros(Xsize, inputDim)];
    Y = zeros(Xsize, 1);
    Y(1) = objFunc(inputs);
    
    for j = 1:replications
        start = 2*(j-1)*inputDim;
        for i = 1:inputDim
            %Compute derivative step
            step = max(cStep*inputs(i), cStep);

            %Compute Positive Step Value
            thisInput = inputs;
            thisInput(i) = thisInput(i) + j*step;
            X(start + 2*i, i+1) = j*step;
            Y(start + 2*i) = objFunc(thisInput);

            %Compute Negative Step Value
            thisInput = inputs;
            thisInput(i) = thisInput(i) - j*step;
            X(start + 2*i + 1, i+1) = -j*step;
            Y(start + 2*i+1) = objFunc(thisInput);
        end
        
    end
    
    %Compute Model and Uncertainty
    %XXinv = inv(X'*X);
    valVect = (X'*X)\X'*Y;
    epsilon = Y - X*valVect;
    eVar = min(var(epsilon), 1E-8);
    cov = eVar*inv(X'*X);
    
%     plot3(X(:,2), X(:,3), Y, '.')
    
end