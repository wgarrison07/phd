
function [minX, minFxEst, minFxUnc] = BoundLineSearch(func, lower, upper, print)
%Initialize Variables
iter = 0;
phi = (sqrt(5) + 1)/2;
range = upper - lower;
fitTol = 1E-2;
sepTol = 1E-4*(upper - lower);
if nargin < 4
    print = false;
end

%Create Input Set
X = zeros(1, 4);
X(1) = lower;
X(2) = upper - range/phi;
X(3) = lower + range/phi;
X(4) = upper;

%Create Function Set
F = zeros(1, 4);
for i = 1:4
    F(i) = func(X(i));
end
avg = mean(abs(F));
count = 4;

%Create Quadratic Fit And Compute Error
[A, err] = FitPoly(F, X);
normErr = err/avg;

if print
    [estMinX, estMinFx] = PolyMin(A, X); 
    fprintf(1, 'Iter=%f, EstMinX=%f, EstMinFx=%f, RelErr=%e\n', iter, estMinX, estMinFx, normErr);
end

%Loop Until Quadratic Approx Meets Error
while normErr > fitTol && (X(4) - X(1)) > sepTol
    iter = iter + 1;
    
    %Update Bounds using Golden Search
    if F(2) < F(3)
        X(4) = X(3);
        F(4) = F(3);
    else
        X(1) = X(2);
        F(1) = F(2);
    end    
    range = X(4) - X(1);
    X(2) = X(4) - range/phi;
    F(2) = func(X(2));
    X(3) = X(1) + range/phi;
    F(3) = func(X(3));
    
    %Update Function Average
    avg = (count*avg + abs(F(2)) + abs(F(3)))/(count+2);
    count = count + 2;
    
    %Recompute Quadratic Approx and Compute Error
    [A, err] = FitPoly(F, X);
    normErr = err/avg;
    
    if print
       [estMinX, estMinFx] = PolyMin(A, X); 
       fprintf(1, 'Iter=%f, EstMinX=%f, EstMinFx=%f, RelErr=%e\n', iter, estMinX, estMinFx, normErr);
    end
end

%Quadratic Curve Meets tolerance, Return Minimum
[minX, minFxEst] = PolyMin(A, X); 
minFxUnc = normErr*avg;

end

function [A, err] = FitPoly(F, X)
%Normalize Input
N = (X - X(1))/(X(4) - X(1));
A = zeros(4,3);

%Estimate F(1) using F(2), F(3), F(4)
A(1,:) = ComputeCoef(2, 3, 4, F, N); 

%Estimate F(2) using F(1), F(3), F(4)
A(2,:) = ComputeCoef(1, 3, 4, F, N); 

%Estimate F(3) using F(1), F(2), F(4)
A(3,:) = ComputeCoef(1, 2, 4, F, N); 

%Estimate F(4) using F(2), F(3), F(4)
A(4,:) = ComputeCoef(1, 2, 3, F, N); 

%Compute Fit Error
est = zeros(1,4);
for i = 1:4
    est(i) = A(i,1) + A(i,2)*N(i) + A(i,3)*N(i)^2;
end
err = abs(est - F);
err = (1/6)*(err(1) + 2*err(2) + 2*err(3) + err(4));

end

function coef = ComputeCoef(i1, i2, i3, F, N)

coef(3) = ((F(i3)-F(i1))/(N(i3)-N(i1)) - (F(i2)-F(i1))/(N(i2)-N(i1)))/(N(i3)-N(i2));
coef(2) = (F(i2)-F(i1))/(N(i2)-N(i1)) - coef(3)*(N(i1)+N(i2));
coef(1) = F(i1) - coef(2)*N(i1) - coef(3)*N(i1)^2;

if abs(coef(3)) < 1E-8
    coef(3) = 1E-8;
end

end

function [estMinX, estMinFx] = PolyMin(A, X)
rng = X(4) - X(1);
lwr = X(1);

estMinXi  = zeros(4,1);
estMinFxi = zeros(4,1);
for i = 1:4
    estMinN      = min(max(-0.5*A(i,2)/abs(A(i,3)),0),1);
    estMinXi(i)  = estMinN*rng + lwr;
    estMinFxi(i) = A(i,1) + A(i,2)*estMinN + A(i,3)*estMinN^2;
end

estMinX  = (1/6)*(estMinXi(1) + 2*estMinXi(2) + 2*estMinXi(3) + estMinXi(4));
estMinFx = (1/6)*(estMinFxi(1) + 2*estMinFxi(2) + 2*estMinFxi(3) + estMinFxi(4));

end