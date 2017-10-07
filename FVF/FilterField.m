%% Vector Field Filtering Function
function field = FilterField(field, input, output, outputCov)

    x = reshape(field.nodeVal',[],1);    
    P = field.covar;
    [valEst, H, valCov] = ComputeEstimateFromField(field, input);
    y = output - valEst;
    S = valCov + outputCov;
    K = P*H'/S;
    x = x + K*y;
    I = eye(length(x));
    ImKH = I - K*H;
    P = ImKH*P; %*ImKH' + K*outputCov*K';
    
    field.nodeVal = reshape(x, field.numNodes, []);
    field.covar = P;
    
end


