close all;

% testH = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2; % + 1E-8*randn(1);
testH = @(x) -(x(1)^2);

varRanges = [-1;1];

field = CreateVectorField(varRanges, 6);

for i = 1:field.numNodes
   startIndx = (i-1)*field.valDim+1;
   indx      = startIndx:startIndx+field.valDim-1;
   [field.nodeVal(i,:), field.covar(indx,indx)] = EvalObjFunc(testH, field.nodePos(i,:), 3); 
end

%Compute Local Error
% 
%Compute Projection Error
x = varRanges(1,1):0.1:varRanges(2,1);
% y = varRanges(1,2):0.1:varRanges(2,2);
% err = zeros(length(x)*length(y),1);
% for i = 1:length(x)
%     for j = 1:length(y)
%         input = [x(i), y(j)];
%         [val, H, ~] = ComputeEstimateFromField(field, input);
%         cov = H*field.covar*H';
%         err((i-1)*length(y) + j) = (val(1) - testH(input))/cov(1,1);
%     end
% end
err = zeros(length(x),1);
for i = 1:length(x)
    input = x(i);
    [val, H, ~] = ComputeEstimateFromField(field, input);
    cov = H*field.covar*H';
    err(i) = (val(1) - testH(input))/cov(1,1);
end

figure();
hist(err);

figure();
hold on;
x = varRanges(1,1):0.01:varRanges(2,1);
plot(field.nodePos(:,1),zeros(length(field.nodePos), 1),'rd');
plot(x, -x.^2);
for i = 1:length(x)
    val = ComputeEstimateFromField(field, x(i));
    plot(x(i), val(1), '.g');
end
% y = varRanges(1,2):0.1:varRanges(2,2);
% [X, Y] = meshgrid(x, y);
% Z = (1 - X.^2) + 100*(Y - X.^2).^2;
% Z = -X.^2 - Y.^2;
% gcf();
% hold off;
% mesh(X, Y, Z);
% hold on;
% plot3(field.nodePos(:,1), field.nodePos(:,2), field.nodeVal(:,1), '+k');
% for i = 1:length(x)
%     for j = 1:length(y)
%         [val, ~, cov] = ComputeEstimateFromField(field, [x(i), y(j)]);
%         plot3(x(i), y(j), val(1) , '.g');
% %         plot3(x(i), y(j), val(1) + sqrt(cov(1,1)), '.b');
% %         plot3(x(i), y(j), val(1) - sqrt(cov(1,1)), '.b');
%     end
% end
% % 
