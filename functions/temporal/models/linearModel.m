function result = linearModel(param, stim, t)
%
% function normrsp = dn_DNmodel(param, stim, t)
% INPUTS  -----------------------------------------------------------------
% params : 2 fields.
%          shift -- time between stimulus onset and when the signal reaches
%          the cortex, in unit of second
%          scale -- response gain.
% OUTPUTS -----------------------------------------------------------------
% normrsp: shifted and scaled linear response.

%% PRE-DEFINED /EXTRACTED VARIABLES

x       = []; % a struct of model parameteres

t_lth   = length(t);
dt      = t(2) - t(1);

normSum = @(x) x./sum(x);


%% SET UP THE MODEL PARAMETERS

fields = {'shift', 'scale'};
x      = toSetField(x, fields, param);


%% COMPUTE THE NORMALIZATION RESPONSE

normrsp = zeros(size(stim, 1),size(stim, 2));
for istim = 1 : size(stim, 1)
    if x.shift > 0
        % ADD SHIFT TO THE STIMULUS -------------------------------------------
        sft       = round(x.shift * (1/dt));
%         stimtmp = [zeros(1,sft) stim(istim, :)]; % hack for now
        stimtmp   = padarray(stim(istim, :), [0, sft], 0, 'pre');
        stim(istim, :) = stimtmp(1 : size(stim, 2));
    end    
    % COMPUTE HTE NORMALIZATION RESPONSE
    normrsp(istim, :) = x.scale.*stim(istim, :);
%     if ismember(istim, round((1:10)/10* size(stim, 1)-1)) % every 10% draw a dot
%         fprintf('[%s]: (%2.0f%%/irf)\n',mfilename,100*(istim/size(stim, 1)));drawnow;
%     end
end
result{1} = normrsp';







end