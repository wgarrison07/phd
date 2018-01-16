function [measEst, measEstCov] = KrigingMeasFunc(PhiInv, Node, Z, ZCov, N, S, R, X)
% PhiInv is an [c, c] matrix that is the inverse of the covariance between
%   the Kriging node points.
% Node is the [c, n] matrix of Kriging node locations
% Z is the [c, 1] vector objective function values at node locations
% ZCov is the [c, c] matrix of covariance of the node Z values
% N is the nugget of the variogram model
% S is the sill (maximum variance) of the variogram model
% R is the range (distance of maximum variance) of the variogram model
% X is an [1, n] vector of the point to be queried

%Compute the Covariance of the query point with node points
c = size(Node, 1);
P = zeros(c, 1);
for i = 1:c
    r = norm(Z(i,:) - X);
    P(i) = (S-N)*exp(-3*r*r/(R*R)) + N;
end

%Compute the Kriging weight
w = PhiInv*P;

%Compute Estimated Objective Function Value
measEst = w*Z;

%Compute Estimated Objective Function Value Uncertainty
measEstCov = w*ZCov*w';

end