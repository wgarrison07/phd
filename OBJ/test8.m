

objFuns = ObjectiveFunctionStructure();

numTrials = 100;
close all;

for i = 1:length(objFuns)
   fun = objFuns{i};
   sinuosity = zeros(numTrials, 1);
   for j = 1:numTrials
        sinuosity(j) = ComputeSinuosity(fun.Func, fun.MinVals, fun.MaxVals); 
   end
%    figure();
%    plot(log(sinuosity));
   fprintf(1,'sinuosity of %s = %f\n', fun.Name, median(sinuosity));
end


% fun = objFuns{1};
% x = 0:0.0001:1;
% y = zeros(size(x));
% for i=1:length(x)
%    y(i) = fun.Func(x(i)); 
% end
% plot(x, y);
% xlabel('x');
% ylabel('f(x)');
% title('Forrester et al. (2008) Function');


% close all;
% objFuns = ObjectiveFunctionStructure();
% obj = objFuns{7};
% step = 0.1;
% x = obj.MinVals(1):step:obj.MaxVals(1);
% y = obj.MinVals(2):step:obj.MaxVals(2);
% [X, Y] = meshgrid(x, y);
% Z = zeros(size(X));
% for i = 1:size(X,1)
%     for j = 1:size(X,2)
%         Z(i, j) = obj.Func([X(i,j), Y(i,j)]);
%     end
% end
% signal = std(Z(:));
% noise = signal / 100;
% noisedZ = Z + randn(size(Z))*noise;
% surf(X, Y, log(Z));
% figure();
% surf(X, Y, noisedZ);