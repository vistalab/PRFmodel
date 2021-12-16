function nidxs = find_stim_idxs(stim, N)
% Helper function for find indices of stimulus step functions.
% 
% INPUTS
%   1) stim: binary stimulus vector of zeros and ones
%   2) N: number of indices to find and return
% 
% OUTPUT
%   idxs: indices of first N stimulus frames in binary vector
% 
% AS 2/2017

if ~isempty(stim)
    idxs = find(stim == 1);
    if nargin < 2
        N = sum(stim == 1);
    end
    N = min([sum(stim == 1) N]);
    nidxs = idxs(1:N);
else
    nidxs = [];
end

end
