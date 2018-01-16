function model = TestDriver()
% rng(2194851);
close all;
obj = ObjectiveFunctionStructure();
f = obj{9};
min = f.MinVals;
max = f.MaxVals;
num = 1000;
x = RandomLatinHypercube([min; max], num);
y = zeros(num,1);
% f.Func = @(x) x.^2;
for i = 1:num
    y(i) = f.Func(x(i,:));
end
sig = var(y);
SNR = 1E+10;
noise = sig/SNR;
fprintf(1, 'Noise Var is: %f\n', noise);
fprintf(1, 'Sig Var is: %f\n', sig);
fprintf(1, 'Truth Nugget: %f\n', noise/(sig + noise));
yn = y + sqrt(noise)*randn(size(y));

model = Kriging(x, yn, 'plot' );

x = RandomLatinHypercube([min; max], num);
y = zeros(num,1);
% f.Func = @(x) x.^2;
for i = 1:num
    y(i) = f.Func(x(i,:));
end

yEst = EvalModel(model, x);
normErr = sqrt(mean((yEst - y).^2))/sqrt(sig);
fprintf(1, 'MVE = %f\n', normErr);

yBar = mean(y);
SStot = sum((y - yBar).^2);
err = yEst - y;
SSres = sum(err.*err);
MRSE = sqrt(SSres/length(err));
R2 = 1 - SSres/SStot;
fprintf(1, 'MRSE = %f\nR2 = %f\n', MRSE, R2);

% EvalModel(model, x);
end