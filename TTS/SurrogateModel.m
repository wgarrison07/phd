
function model = SurrogateModel()

%Basic Values
model = struct();
model.Type = '';
model.InputDim = 0;
model.OutputDim = 1;
model.InputMu = 0;
model.InputSigma = 1;
model.OutputMu = 0;
model.OutputSigma = 1;

%Fit Staticstics
model.RMSE = 0;
model.R2 = 0;
model.AdjR2 = 0;
model.MVE = 0;

end



