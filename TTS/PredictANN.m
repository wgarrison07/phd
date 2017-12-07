function [yEst, X] = PredictANN(W1, W2, x)
    n = size(x, 1);
    yEst = zeros(n, size(W2, 1));
    X = zeros(n, size(W2, 2));
    x = [x, ones(n,1)];
    for i = 1:n
        z = W1*x(i, :)';
        g = [tanh(z);1];
        X(i, :) = g';
        yEst(i, :) = (W2*g)';
    end
end