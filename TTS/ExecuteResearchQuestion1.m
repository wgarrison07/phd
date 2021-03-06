function ExecuteResearchQuestion1()

% types = {'PolyReg', 'SplineReg', 'RBF', 'Kriging', 'ANN'};
% types = {'PolyReg', 'SplineReg', 'ANN'};
SNRs  = [1E8, 100, 10, 5, 1, 0.1];
OBJs  = ObjectiveFunctionStructure()';

%Create Directory for Data to Save
dirName = 'ResearchQuestionOne';
mkdir(fullfile('Results',dirName));


i = 0;
n = length(types)*length(SNRs)*length(OBJs);
timer = tic();
for iOBJ = OBJs
    for iType = types
        for iSNR = SNRs
            %Test Model and Save Results
           [R2, samp] = ExecuteModelTrainingTest(iOBJ{1}, iType{1}, iSNR, i, n, timer); 
           fileName = sprintf('%s_%s_%4.1f.mat', iOBJ{1}.Name, iType{1}, iSNR);
           save(fullfile('Results',dirName, fileName), 'R2', 'samp');
           i = i+1;
        end    
    end
end

end