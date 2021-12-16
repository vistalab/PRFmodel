
function output = convCut2(stimulus, impulse, nTerms)
%
% INPUTS -----------------------------------------------------
% sconvCut(stimulus, impulse, nTerms)
%
% INPUTS -----------------------------------------------------
% stimulus : a stimulus time course
% impulse  : an impulse response function
% nTerms   : number of terms after cutting
%
% OUTPUT(S) --------------------------------------------------
% output   : output cutted between 1 and nTerms

% % DEPENDENCIES ----------------------------------------------

%%

output = conv(stimulus, impulse);

output = output(1:nTerms);


end