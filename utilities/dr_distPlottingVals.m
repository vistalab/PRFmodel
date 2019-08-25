function [y, mu, sigma, mn, mx] = dr_distPlottingVals(x)
    mu    = mean(x, 'omitnan');
    sigma = std(x, 'omitnan');
    mn    = mu - 5 * sigma;
    mx    = mu + 5 * sigma;
    y     = mn:((mx-mn)/1000):mx;
end