function output = trf_simpNormModel(params, stimulus, t)

% INPUTS:
% params : 
%       sigma
%       n
%       tau1

% OUTPUTS:

%% default parameters

n = 2;

%% useful functions

normMax = @(x) x./max(x);

%% compute impulse response

irf    = gammaPDF(t, params.tau1, 2);

%% compute normalization response

for k = size(stimulus, 1)
    linRsp(k, :) = convCut(irf, stimulus(k, :), length(stimulus));
    num(k, :)    = linRsp.^n;
    den(k, :)    = params.sigma.^n + num(k, :);
    output(k, :) = params.scale.*num(k, :)./den(k, :);
end

end