%Compute Validation Error of a Model
function [err, C] = VE(X, y, Xval, yval, lambda)

%Update Lambda to ensure positive definiteness
I = eye(size(X,2));
L = max(lambda, 1E-3)*I;
XtX = X'*X;
[cholX, p] = chol(XtX + L);
while p ~= 0
    L = 2*L;
    [cholX, p] = chol(XtX + L);
end  

temp = cholX\I;
C = (temp*temp')*(X*y);
yEst = Xval*C;
err = norm(yEst - yval);

end