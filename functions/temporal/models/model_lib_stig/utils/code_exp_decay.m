function resp_decay = code_exp_decay(resp_in, starts, stops, decay_exp, fs)
% Helper function for coding stimulus-specific exponential response decay. 
%
% INPUTS
%   1) resp_in: input activity matrix (frames x predictors)
%   2) starts: beginnings of decay activity windows (seconds)
%   3) stops: ends of decay activity windows (seconds)
%   4) decay_exp: exponential function modeling decay of activity
%   5) fs: temporal sampling rate of resp_in and resp_out (Hz)
%
% OUTPUT
%   resp_decay: output activity matrix with decay (frames x predictors)
%
% AS 10/2017

if ~isempty(resp_in)
    decay_fun = zeros(size(resp_in, 1), 1);
    for ss = 1:length(starts)
        start_idx = round(starts(ss) * fs); stop_idx = round(stops(ss) * fs);
        decay_idxs = start_idx + 1:stop_idx; dl = length(decay_idxs);
        decay_idxs = decay_idxs(1:min([dl length(decay_exp)]));
        decay_fun(decay_idxs) = decay_exp(1:min([dl length(decay_exp)]));
    end
    resp_decay = resp_in .* repmat(decay_fun, 1, size(resp_in, 2));
else
    resp_decay = [];
end

end

