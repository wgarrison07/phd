%Generate Overfitting Example
close all;
rng(2194851);
noise = 3;

%Generate Forrester Data
fx = (0:0.06:1)';
fy = zeros(size(fx));
ny = zeros(size(fx));
for i = 1:length(fx)
    fy(i) = Forrester(fx(i));
    ny(i) = fy(i) + noise*randn();
end
snr = var(fy)/noise^2;
fprintf(1, 'SNR = %f\n', snr);

n = length(fx);
x = (0:0.005:1)';
X = zeros(n, n);
X2 = zeros(length(x), n);
for i = 1:n
    X(:, i) = fx.^(i-1);
    X2(:, i) = x.^(i-1);
end

C = (X'*X)\(X'*ny);
yEst = X2*C;

figure();
plot(fx, fy, '-g');
hold on;
plot(fx, ny, '.k');
plot(x, yEst, '-r');
ylim([-10, 25]);


model = PolyReg(fx, ny);
yEst = EvalModel(model, x);
figure();
plot(fx, fy, '-g');
hold on;
plot(fx, ny, '.k');
plot(x, yEst, '-r');