function pm = timeSeriesCompute(pm)
% Compute the RF from the current parameters
%
%
% See also


% BOLD signal
% Predicted prf response
tSamples = size(pm.stimulus.binary,3);
pm.BOLD.timeSeries = zeros(1,tSamples);
for tt = 1:tSamples
    % This is called the hadamard product.  It is the pointwise
    % multiplication of the RF with the stimulus.  The hadProduct is
    % the same size as the stimulus
    hadProduct = pm.stimulus.binary(:,:,tt) .* pm.RF.values;
    
    % Now, we add up all of the hadProduct values
    pm.BOLD.timeSeries(tt) = sum(hadProduct(:));
end

% pm.timeSeries is the signal prior to convolution

% Convolution between the timeSeries and the HRF
pm.BOLD.predicted = conv(pm.BOLD.timeSeries, ...
                         pm.HRF.values, ...
                         'same');


end