close all;
% rng(319954);
%Generate Forrester Data
fx = 0:0.01:1;
fy = zeros(size(fx));
for i = 1:length(fx)
    fy(i) = Forrester(fx(i));
end

figure();
subplot(1, 2, 1);
noise = 10;
x = [0.1, 0.3, 0.5, 0.7, 0.9];
zx = zeros(5*length(x), 1);
zy = zeros(5*length(x), 1);
count = 1;
y = zeros(size(x));
for i = 1:length(x)
    for j = 1:10
        zx(count) = x(i);
        zy(count) = Forrester(x(i)) + noise*randn();
        y(i) = y(i) + zy(count);
        count = count + 1;
    end
    y(i) = y(i)/10;
end
plot(zx, zy, '.k');
hold on;
plot(x, y, 'b+');
plot(fx, fy, '-r');
X = [ones(5, 1), x', x'.^2, x'.^3, x'.^4];
C = (X'*X)\(X'*y');
X = [ones(size(fx')), fx', fx'.^2, fx'.^3, fx'.^4];
my = X*C;
% model = PolyReg(x', y', 'order', 5);
% my = EvalModel(model, fx');
plot(fx, my, '-g');
err = sqrt(mean((fy' - my).^2));
title(sprintf('Root Mean Squared Error = %f', err));

subplot(1, 2, 2);
x = 0:(1/49):1;
y = zeros(size(x));
for i = 1:length(x);
    y(i) = Forrester(x(i)) + noise*randn();
end
plot(x, y, '.k');
hold on;
plot(fx, fy, '-r');
model = PolyReg(x', y');
my = EvalModel(model, fx');
plot(fx, my, '-g');
err = sqrt(mean((fy' - my).^2));
title(sprintf('Root Mean Squared Error = %f', err));
