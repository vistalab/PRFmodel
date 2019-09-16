function [y, mu, sigma, mn, mx, med, ci] = dr_distPlottingVals(x)
    mu    = mean(x, 'omitnan');
    sigma = std(x, 'omitnan');
    mn    = mu - 5 * sigma;
    mx    = mu + 5 * sigma;
    y     = mn:((mx-mn)/1000):mx;
    
    
    % Obtain the median and the 25%-75% confidence intervals
    med    = median(x, 'omitnan');
    % Define the required confidence intervals
    CIrange        = 50;
    twoTailedRange = (100 - CIrange)/2;
    % ci             = prctile(x, [twoTailedRange, 100-twoTailedRange]);
    ci             = quantile(x, [twoTailedRange/100, (100-twoTailedRange)/100]);
    
    
end