function model = TestDriver()
close all;
obj = ObjectiveFunctionStructure();
f = obj{3};
min = f.MinVals;
max = f.MaxVals;
num = 200;
x = RandomLatinHypercube([min; max], num);
y = zeros(num,1);
for i = 1:num
    y(i) = f.Func(x(i,:));
end
sig = var(y);
SNR = 20;
noise = sig/SNR;
y = y + sqrt(noise)*randn(size(y));


% y2 = zeros(num,1);
% for i = 1:num
%     y2(i) = f.Func(x(i,:)) + 4*randn();
% end
% 
% 
% y3 = zeros(num,1);
% for i = 1:num
%     y3(i) = f.Func(x(i,:)) + 4*randn();
% end
% model = RBF([x;x;x], [y;y2;y3]);


model = Kriging(x, y);
end