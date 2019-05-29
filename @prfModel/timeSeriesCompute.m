function pm = timeSeriesCompute(pm)
% Compute the RF from the current parameters
%
%
% See also


% BOLD signal
% Predicted prf response

% Obtain the values if it is a path
stimValues     = stimValuesRead(pm.stimulus.values);
% Obtain number of samples    
tSamples = size(stimValues,3);

% Initialize timeSeries
pm.BOLD.timeSeries = zeros(1,tSamples);
for tt = 1:tSamples
    % This is called the hadamard product.  It is the pointwise
    % multiplication of the RF with the stimulus.  The hadProduct is
    % the same size as the stimulus
    hadProduct = stimValues(:,:,tt) .* pm.RF.values;
    
    % Now, we add up all of the hadProduct values
    pm.BOLD.timeSeries(tt) = sum(hadProduct(:));
end

% pm.timeSeries is the signal prior to convolution

% Convolution between the timeSeries and the HRF
pm.BOLD.predicted = conv(pm.BOLD.timeSeries, ...
                         pm.HRF.values, ...
                         'same');
pm.BOLD.tSamples  = 1:pm.TR:pm.TR*tSamples;


end