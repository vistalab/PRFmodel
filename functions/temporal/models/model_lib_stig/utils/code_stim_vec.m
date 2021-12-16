function stim_out = code_stim_vec(stim_in, idxs, cond, val)
% Helper function for coding stimulus step functions.
% 
% INPUTS
%   1) stim_in: stimulus design matrix (time x condition)
%   2) idxs: row indicdes to code
%   3) cond: column index to code
%   4) val: value to code at idxs in cond
% 
% OUTPUT
%   stim_out: stimulus design with idxs in cond set to val
% 
% AS 2/2017

if nargin < 3; cond = 1; end
if nargin < 4; val = 1; end
stim_out = [];
if ~isempty(stim_in)
    stim_out = stim_in;
    stim_out(idxs, cond) = val;
end

end
