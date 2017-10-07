%% Vector Field Creation Function
function vectorField = CreateVectorField(varRanges, numNodes)

vectorField = struct();
inputDim = size(varRanges,2);
vectorField.inputRange = varRanges(2,:)-varRanges(1,:);
vectorField.inputDim = inputDim;
vectorField.valDim = inputDim + 1;
vectorField.numNodes = numNodes;
vectorField.nodePos = RandomLatinHypercube(varRanges, numNodes);
vectorField.nodeVal = zeros(numNodes, (inputDim+1));
vectorField.covar   = 1E10*eye(numNodes*(inputDim+1));
%vectorField.nodeIndx = delaunayn(vectorField.nodePos);


end