
function output = convCut(stimulus, impulse, nTerms)
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

output = conv(squeeze(stimulus), squeeze(impulse), 'full');

output = output(1:nTerms);


end