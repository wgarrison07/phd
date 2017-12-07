function PlotModel(model)

x = (-3:0.1:3).*model.InputSigma(1) + model.InputMu(1);
y = (-3:0.1:3).*model.InputSigma(1) + model.InputMu(1);
[X, Y] = meshgrid(x, y);
Z = zeros(size(X));
for i = 1:size(Z, 1)
    for j = 1:size(Z, 2)
        Z(i, j) = EvalModel(model, [X(i, j), Y(i, j)]);
    end
end

Z = (Z - model.OutputMu) / model.OutputSigma;
X = (X - model.InputMu(1)) / model.InputSigma(1);
Y = (Y - model.InputMu(2)) / model.InputSigma(2);
surf(X, Y, Z);

end