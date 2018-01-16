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
        hold on;
    end
        
    %Load File and Plot Data
    load(fullfile(folder, files(i).name), 'R2', 'samp');
    plot(samp, R2, '-b');
    title(plotname);
    
%     C = [(1-R2(1)); 0];
%     t = samp';
%     err = C(1)./(C(2) + t) - (1-R2');
%     J = zeros(length(samp), 2);
%     for j = 1:10
%         J(:, 1) = 1./(C(2) + t);
%         J(:, 2) = -C(1)./(C(2) + t).^2;
%         C = C + ((J'*J)\(J'*err))
%         err = C(1)./(C(2) + t) - (1-R2');
%     end
        
    X = [ones(size(samp')), 1./(samp').^0.5];
    C = (X'*X)\(X'*(1-R2'));
    R2est = 1 - X*C;
    plot(samp, R2est, '--g');
    
%     axis([10, 500, -0.1, 1]);

end
end