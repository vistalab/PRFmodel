function cpred = convolve_vecs_equal(pred, irf)
% Convolve predictors with HRF, resample to temporal resolution of fMRI
% measurements, and clip extra frames extending beyond measurement period.
% Note that this implementation is much faster than using conv2.
% 
% INPUTS
%   1) pred: predictors to be convolved with impulse response function
%   2) irf: impulse response function (sampling rate matched to pred)
% 
% OUTPUT
%   cpred: convolved predictors without resample
% 
% this is without resampling

[nframes, npreds] = size(pred);
cpred = zeros(ceil(nframes), npreds);
for pp = 1:npreds
    cs = fconv(irf, pred(:, pp));
    cpred(:, pp) = cs(1:nframes);
end

end

