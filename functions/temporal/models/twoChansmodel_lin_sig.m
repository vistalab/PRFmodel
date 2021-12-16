function result = twoChansmodel_lin_sig(param, stimulus, t)
% function output = dn_2Chansmodel(params, stimulus, t, dt)

% INPUTS : 
% params have fields for two beta weights, b1 and b2

%% useful functions

% makeIRF = @(A, B, C, t)(t/A).^8 .* exp(-t/A) - 1 / B .* (t/C).^9 .* exp(-t/C);
dt      = t(2) - t(1);
normSum = @(x) x./sum(x);
normMax = @(x) x./max(x);
normRange = @(x) (x-min(x))/(max(x)-min(x));

%% SET UP THE MODEL PARAMETERS
x       = [];
fields = {'tau_s', 'tau_ae', 'lambda_p', 'kappa_p', 'kappa_n', 'weight','shift'};
x      = toSetField(x, fields, param);

% param = [4.93 10000 0.1 0.3 3 0.5 0];

% default for exponential time constants & decay
% tau_s = 4.93;
% tau_ae = 10000;

% default for sigmoid functions
% lambda_p = .1; kappa_p = 3; kappa_n = 3;

%% defualt parameters


if nargin ~= 3; error('Unexpected input arguements.'); end

% default paramters for channel IRFs
dt      = t(2) - t(1);

n1 = 9; n2 = 10; kappa = 1.33; fs=1/dt;

% sustained and transient weight
% weight = 0.5;
%% make impulse response function 

nrfS = tch_irfs('S', x.tau_s, n1, n2, kappa, fs);
nrfT = tch_irfs('T', x.tau_s, n1, n2, kappa, fs);

% nrfS = normMax(normSum(nrfS));
% nrfT = normMax(normSum(nrfT));
% t_irf   = 1 : dt : 60000;

adapt_exp = exp(-(1:60000) / x.tau_ae);

%% make impulse response function

adapt_acts = zeros(size(stimulus, 1),size(stimulus, 2));
predTs = adapt_acts; output=adapt_acts;
for k = 1 : size(stimulus, 1)
    
    % ADD SHIFT TO THE STIMULUS -------------------------------------------
    if x.shift > 0
        sft       = round(x.shift * (1/dt));
        stimtmp   = padarray(stimulus(k, :), [0, sft], 0, 'pre');
        stimulus(k, :) = stimtmp(1 : size(stimulus, 2));
    end
    
    % code stimulus on and offs -------------------------------------------
%     [onsets,offsets, ~] = st_codestim(stimulus(k,:)',fs);

    % sustained with adaptation -------------------------------------------
    adapts = convCut(nrfS, stimulus(k, :), length(stimulus));
%     adapt_acts(k,:) = code_exp_decay(adapts', onsets,offsets,adapt_exp,fs);
    adapt_acts(k,:) = x.weight * adapts(k,:);
    adapt_acts(k,:) = normMax(adapt_acts(k,:));

    % transient with sigmoid -------------------------------------------
    predT = convCut(nrfT, stimulus(k, :), length(stimulus));
    predTs(k, :) = tch_sigmoid(predT, x.lambda_p, x.kappa_p, x.lambda_p, x.kappa_n);
    predTs(k, :) = (1-x.weight) * predTs(k,:);
    predTs(k, :) = normMax(predTs(k, :));
% 


    % weighted sum of the two chennels -------------------------------------------
%     scaleFactor = max(adapt_acts(k,:))/max(predTs(k,:));
%     predTs(k,:) = scaleFactor * predTs(k,:);
%     output(k, :) = x.weight * adapt_acts(k,:) +  (1-x.weight) * predTs(k,:);

%     output(k, :) = adapt_acts(k,:) + predTs(k,:);

%     if ismember(k, round((1:10)/10* size(stimulus, 1)-1)) % every 10% draw a dot
%         fprintf(1,'[%s]: (%2.0f%%/irf)\n',mfilename,100*(k/size(stimulus, 1)));drawnow;
%     end
end

% sustained, transient, resp (weighted sum)
result{1}  = adapt_acts';
result{2}  = predTs';
% result{3}  = output';





end