function predictedTS = pmNoiseWhite(predictedTS, k)
%   Adds white noise to a time series
    n = k * mean(predictedTS);
    noise = n * randn(size(predictedTS));
    predictedTS = predictedTS + noise;
end

