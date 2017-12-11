function GenerateResearchQuestionOnePlots()

%Get Contents of Directory
folder = fullfile('Results', 'ResearchQuestionOne');
files = dir(folder);

%Create Figure Mapping Structure and Plot Files
close all;
plots = struct();
set(0,'DefaultFigureWindowStyle','docked')
for i = 3:length(files)
    %Parse File
    [~, name, ext] = fileparts(files(i).name);
    if ~strcmp(ext, '.mat')
        continue;
    end
    nameParts = strsplit(name, '_');
    obj = nameParts{1};
    type = nameParts{2};
    snr  = str2num(nameParts{3});
    plotname = [obj, type];
    
    %Get Active Figure
    if ismember(plotname, fieldnames(plots))
        figure(plots.(plotname).h);
        hold on;
    else
        h = figure();
        plots.(plotname).h = h;
    end
        
    %Load File and Plot Data
    load(fullfile(folder, files(i).name), 'MVE', 'samp');
    loglog(samp, MVE);
    title(plotname);
    axis([10, 500, 1E-3, 3E1]);

end
end