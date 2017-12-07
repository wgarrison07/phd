function var = VarianceModel(N, S, lag, Type)
    lag = max(min(lag, 1), 0);
    if strcmp(Type, 'linear')
        var = (S-N)*lag + N;
    elseif  strcmp(Type, 'spherical')
        var = (S-N)*(3/2*lag - 0.5*lag.^3) + N;
    elseif strcmp(Type, 'exponential')
        var = (S-N)*(1 - exp(-3*lag)) + N;
    elseif strcmp(Type, 'gaussian')
        var = (S-N)*(1 - exp(-3*lag.^2)) + N;
    else
        error('Unsupported Model Type');
    end
end