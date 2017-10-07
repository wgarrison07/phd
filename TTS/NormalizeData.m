

function [data, mu, sigma] = NormalizeData(data, mu, sigma)

if nargin < 2
    mu = mean(data, 1);
    sigma = std(data, 1);
end

n = size(data,1);
data = (data - repmat(mu,n,1))./repmat(sigma, n, 1);

end