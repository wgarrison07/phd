function [valEst, H, valCov] = ComputeEstimateFromField(field, input)    
    %Compute Distance for all nodes from input for scaling
    nodeDist = repmat(input, field.numNodes, 1) - field.nodePos;
    distStd = std(nodeDist);
    normDist = norm(field.inputRange);
    
    
    %Compute the weight for each node (pure distance)
%     nodeWeight = zeros(field.numNodes, 2);
%     for i = 1:field.numNodes
%         nodeWeight(i, 1) = exp(norm(nodeDist(i,:))^-2);
%         nodeWeight(i, 2) = nodeWeight(i, 1);
%     end
    
    %Compute the weight or each node component (uncertainty * distance)
    nodeWeight = zeros(field.numNodes, 2);
    for i = 1:field.numNodes
%         startIndx = (i-1)*field.valDim+1;
%         indx = startIndx:startIndx+field.inputDim;
%         uncert = diag(field.covar(indx, indx));
%         uncertGrad = uncert(2:end)'*abs(nodeDist(i,:))';
%         uncertPos  = uncert(1) + uncertGrad;
        distFactor = norm(nodeDist(i,:))/normDist;
         nodeWeight(i, 1) = min(1E16, exp(-3*norm(nodeDist(i,:))^2/normDist));
        nodeWeight(i, 2) = min(1E16, exp(-distFactor));
    end
    
    hist(nodeWeight(1));
    
    %Normalize Weights
    sumWeight = sum(nodeWeight);
    nodeWeight = nodeWeight./repmat(sumWeight,field.numNodes,1);
    sumWeight = sum(nodeWeight);
    nodeWeight = nodeWeight./repmat(sumWeight,field.numNodes,1);

%     figure();
%     hist(nodeWeight);
    
    %Compute Measurement Transform from Weights
    H = zeros(field.valDim, numel(field.nodeVal));
    nodeEst = zeros(field.numNodes, field.valDim);
    for i = 1:field.numNodes
        %Compute flattening indices
        startIndx = (i-1)*field.valDim+1;
        gradIndx = (startIndx+1):(startIndx+field.inputDim);
        
        %Assign position estimate weights
        H(1, startIndx) = nodeWeight(i,1);
        H(1, gradIndx) = nodeWeight(i,1)*nodeDist(i,:);
        
        %Assign gradient estimate weights
        for j = 1:field.inputDim
            H(j+1, startIndx+j) = nodeWeight(i,2);
        end        
        
        %Compute weighted estimate
        nodeEst(i, 1) = nodeWeight(i,1) * (field.nodeVal(i,1) + field.nodeVal(i,2:end)*nodeDist(i,:)');
        nodeEst(i, 2:end) = nodeWeight(i,2) * field.nodeVal(i,2:end);
    end
    
    %Compute Observation Measurement
    valEst = H*reshape(field.nodeVal',[],1);
    
    %Compute Weighted Covariance of Observation
    valCov = zeros(field.valDim, field.valDim);
    for j = 1:field.valDim
        for k = 1:field.valDim
            sumWeight = 0;
            for i = 1:field.numNodes
                valCov(j,k) = valCov(j,k) + nodeWeight(i,min(j,2))*nodeWeight(i,min(k,2))*(nodeEst(i,j) - valEst(j)) * (nodeEst(i,k) - valEst(k));
                sumWeight = sumWeight + nodeWeight(i,min(j,2))*nodeWeight(i,min(k,2));
            end
            valCov(j,k) = valCov(j,k) / sumWeight;
        end
    end
    
    
end
