
testH = @(x) x(1)^2 + x(2)^2 + 1E-6*randn(1);

% [val, cov] = EvalObjFunc(testH, [0;0]);

FVF(testH, [-1, -1; 1, 1]);
