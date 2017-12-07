function model = TestDriver()
close all;
obj = ObjectiveFunctionStructure();
f = obj{11};
min = f.MinVals;
max = f.MaxVals;
num = 500;
x = RandomLatinHypercube([min; max], num);
y = zeros(num,1);
% f.Func = @(x) x.^2;
for i = 1:num
    y(i) = f.Func(x(i,:));
end
sig = var(y);
SNR = 1E10;
noise = sig/SNR;
fprintf(1, 'Noise Var is: %f\n', noise);
fprintf(1, 'Sig Var is: %f\n', sig);
fprintf(1, 'Truth Nugget: %f\n', noise/(sig + noise));
yn = y + sqrt(noise)*randn(size(y));

model = RBF(x, yn, 'plot');

yEst = EvalModel(model, x);
normErr = sqrt(mean((yEst - y).^2))/sqrt(sig);
fprintf(1, 'MVE = %f\n', normErr);

% EvalModel(model, x);
end