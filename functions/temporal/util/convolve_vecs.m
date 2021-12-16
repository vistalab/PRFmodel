function cpred = convolve_vecs(pred, irf, fs_in, fs_out)
% Convolve predictors with HRF, resample to temporal resolution of fMRI
% measurements, and clip extra frames extending beyond measurement period.
% Note that this implementation is much faster than using conv2.
% 
% INPUTS
%   1) pred: predictors to be convolved with impulse response function
%   2) irf: impulse response function (sampling rate matched to pred)
%   3) fs_in: temporal sampling rate of input predictor (Hz)
%   4) fs_out: temporal sampling rate of output predictor (Hz)
% 
% OUTPUT
%   cpred: convolved predictors (resampled to fs_out)
% 
% AS 2/2017

[nframes, npreds] = size(pred); lmi_in = 1; lmi_out = 1;
cpred = zeros(ceil(nframes / (fs_in * (1 / fs_out))), npreds);
if rem(fs_in, 1) ~= 0; lmi_in = double(lcm(sym(fs_in), 1) / fs_in); end
if rem(fs_out, 1) ~= 0; lmi_out = double(lcm(sym(fs_out), 1) / fs_out); end
P = fs_out * lmi_in * lmi_out; Q = fs_in * lmi_in * lmi_out;
for pp = 1:npreds
    cs = fconv(irf, pred(:, pp));
    cpred(:, pp) = resample(cs(1:nframes), P, Q);
end

end

