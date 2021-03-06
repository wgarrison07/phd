function PlotModel(model)

x = (-1:0.05:1).*model.InputSigma(1) + model.InputMu(1);
y = (-1:0.01:1).*model.InputSigma(1) + model.InputMu(1);
[X, Y] = meshgrid(x, y);
out = EvalModel(model, [X(:), Y(:)]);
Z = reshape(out, size(X));
Z = (Z - model.OutputMu) / model.OutputSigma;
X = (X - model.InputMu(1)) / model.InputSigma(1);
Y = (Y - model.InputMu(2)) / model.InputSigma(2);
surf(X, Y, Z);

end