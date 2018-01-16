close all;
rng(319954);
noise = 3;

%Generate Forrester Data
fx = (0:0.01:1)';
fy = zeros(size(fx));
for i = 1:length(fx)
    fy(i) = Forrester(fx(i));
end
snr = var(fy)/noise^2;
fprintf(1, 'SNR = %f\n', snr);

model1 = PolyReg(fx, fy);
my = EvalModel(model1, fx);
figure();
subplot(1, 2, 1);
plot(fx, fy, '.k');
hold on;
plot(fx, my, '-g');
plot(fx, fy, '-r');
err = sqrt(mean((fy - my).^2));
title(sprintf('Root Mean Squared Error = %f', err));
model1.C

subplot(1, 2, 2);
noisedfy = fy + noise*randn(size(fy));
model2 = PolyReg(fx, noisedfy);
my = EvalModel(model2, fx);
plot(fx, noisedfy, '.k');
hold on;
plot(fx, my, '-g');
plot(fx, fy, '-r');
err = sqrt(mean((fy - my).^2));
title(sprintf('Root Mean Squared Error = %f', err));
model2.C
