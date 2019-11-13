function diff = pm_fitSinusoidal(tSeries, period, varargin)
%PM_FITSINUSOIDAL Fit sinusoidal knowing the period and return the difference


    % Load example data
    % tSeries = realNormtsmn(:,3);
    % period  = 16;
    
    
    % Read the inputs
    varargin = mrvParamFormat(varargin);
    p = inputParser;
    p.addRequired('tSeries', @isnumeric);
    p.addRequired('period' , @isnumeric);
    p.addParameter('plotit', false,@islogical);
    p.parse(tSeries,period,varargin{:});
    
    plotit = p.Results.plotit;
    
    % Check the inputs
    if length(size(tSeries))>2
        error('The input to this function needs to be a vector')
    end
    if size(tSeries,1) < size(tSeries,2)
        tSeries = tSeries';
    end
    if length(period) ~= 1
        error('period needs to be a single number, not a vector')
    end
    
    % Obtain support variables
    nTimes = size(tSeries,1);
    
    % Do the fit
    t      = (1:nTimes)';
    X      = ones(nTimes,3);
    X(:,2) = cos((2*pi)/period*t);
    X(:,3) = sin((2*pi)/period*t);
    y      = tSeries;
    beta   = X\y;
    yhat   = beta(1)+beta(2)*cos((2*pi)/period*t)+beta(3)*sin((2*pi)/period*t);
    if plotit
        mrvNewGraphWin;
        plot(t,y,'b');
        hold on
        plot(t,yhat,'r','linewidth',2);
    end
    % Calculate the difference (noise) and return it
    diff = y - yhat;
end

