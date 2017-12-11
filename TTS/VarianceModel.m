function var = VarianceModel(N, S, R, dist, Type)
    R = max(R, 1E-8);
    if strcmpi(Type, 'linear')
        dist = max(min(dist, R), 0);
        var = (S-N)*dist/R + N;
    elseif  strcmpi(Type, 'spherical')
        dist = max(min(dist, R), 0);
        var = (S-N)*(3*dist/(2*R) - 0.5*dist.^3/R^3) + N;
    elseif strcmpi(Type, 'exponential')
        var = (S-N)*(1 - exp(-3*dist/R)) + N;
    elseif strcmpi(Type, 'gaussian')
        var = (S-N)*(1 - exp(-3*dist.^2/R^2)) + N;
    else
        error('Unsupported Model Type');
    end
end