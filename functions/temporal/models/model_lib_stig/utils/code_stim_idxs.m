function idxs = code_stim_idxs(starts, stops)
% Helper function for coding indices of stimulus step functions.
% 
% INPUTS
%   1) starts: row vector of indices coding the starts of stimulus events
%   2) stops: row vector of indices coding the ends of the events
% 
% OUTPUT
%   idxs = [starts(1):stops(1) starts(2):stops(2) ... starts(n):stops(n)]
% 
% AS 2/2017

lens = stops - starts + 1;
csum = [lens(:) repmat(-1, numel(lens), max(lens) - 1)]';
csum = cumsum(csum) > 0;
strk = double(csum);
strk(1, :) = starts;
strk = cumsum(strk);
idxs = strk(csum);

end
