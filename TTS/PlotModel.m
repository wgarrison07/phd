function PlotModel(model)

x = (-3:0.1:3).*model.InputSigma(1) + model.InputMu(1);
y = (-3:0.1:3).*model.InputSigma(1) + model.InputMu(1);
[X, Y] = meshgrid(x, y);
out = EvalModel(model, [X(:), Y(:)]);
Z = reshape(out, size(X));
Z = (Z - model.OutputMu) / model.OutputSigma;
X = (X - model.InputMu(1)) / model.InputSigma(1);
Y = (Y - model.InputMu(2)) / model.InputSigma(2);
surf(X, Y, Z);

end