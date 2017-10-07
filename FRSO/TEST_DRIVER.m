
close all;

range = [-32, -32; 32, 32];

% testH = @(x) 2*x(1)^2 - 1.05*x(1)^4 + x(1)^6/6 + x(1)*x(2) + x(2)^2;
testH = @(x) -20*exp(-0.2*sqrt(0.5*sum(x(1:2).^2))) - exp(0.5*(cos(2*pi()*x(1)) + cos(2*pi()*x(2)))) + 20 + exp(1) + 0*randn();
% testH = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2 + 1E2*randn(1);
% testH = @(x) x(1)^2 + x(2)^2;

% [x, f, exit, output] = fminsearch(testH, [0; 0]);

% inputDim = 2;
% maxOrder = 2;
% numCoef = (maxOrder + 1)^inputDim;
% coef = zeros(numCoef, 1);
% coefCov = 1E10*eye(numCoef);
% powerMat = zeros(numCoef, inputDim);
% line = 1;
% powers = zeros(1, inputDim);
% while line <= numCoef
%     powerMat(line,:) = powers;
%     for i = 1:inputDim
%         if powers(i) == maxOrder
%             powers(i) = 0;
%         else
%             powers(i) = powers(i) + 1;
%             break;
%         end
%     end
%     line = line + 1;
% end
% 
% 
% coef(3) = 1;
% coef(7) = 1;
% pos = [2;2];
% x = ComputeEstimate(pos, coef, powerMat);
% H = LinearizeObservation(pos, coef, powerMat);
% xp = H*coef;
% 
% 
% blarg = 12;


 FRSO(testH, range);