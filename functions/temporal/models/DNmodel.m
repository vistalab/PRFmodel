function result = DNmodel(param, stim, t)
%
% function normrsp = dn_DNmodel(param, stim, t)
% INPUTS  -----------------------------------------------------------------
% params : 6 fields.
%          tau1 -- irf peak time, in unit of second
%          weight -- the weight in the biphasic irf function, set weight to
%          0 if want to use uniphasic irf function.
%          tau2 -- time window of adaptation, in unit of second
%          n -- exponent
%          sigma -- semi-saturation constant
%          shift -- time between stimulus onset and when the signal reaches
%          the cortex, in unit of second
%          scale -- response gain.
% OUTPUTS -----------------------------------------------------------------
% normrsp: normalization response.

% 04/05 I added parameter field "n", 

%% PRE-DEFINED /EXTRACTED VARIABLES

x       = []; % a struct of model parameteres

t_lth   = length(t);
dt      = t(2) - t(1);

normSum = @(x) x./sum(x);

%% SET UP THE MODEL PARAMETERS

fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
x      = toSetField(x, fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

% HERE I ASSUME THAT THE NEGATIVE PART OF THE IMPULSE RESPONSE HAS A TIME
% CONSTANT 1.5 TIMES THAT OF THE POSITIVE PART OF THE IMPULSE RESPONSE
if x.tau1 > 0.5, warning('tau1>1, the estimation for other parameters may not be accurate'); end
    
t_irf   = dt : dt : 5;

irf_pos = gammaPDF(t_irf, x.tau1, 2);
irf_neg = gammaPDF(t_irf, x.tau1*1.5, 2);
irf     = irf_pos - x.weight.* irf_neg;

%% COMPUTE THE DELAYED REPSONSE FOR THE NORMALIZATION POOL

irf_norm = normSum(exp(-t_irf/x.tau2));
normRange = @(x) (x-min(x))/(max(x)-min(x));

%% COMPUTE THE NORMALIZATION RESPONSE

linrsp = zeros(size(stim, 1),size(stim, 2));
numrsp = linrsp; poolrsp=linrsp; demrsp = linrsp; normrsp = linrsp;
for istim = 1 : size(stim, 1)
    
    if x.shift > 0
        % ADD SHIFT TO THE STIMULUS -------------------------------------------
        sft       = round(x.shift * (1/dt));
        stimtmp   = padarray(stim(istim, :), [0, sft], 0, 'pre');
        stim(istim, :) = stimtmp(1 : size(stim, 2));
    end
    
    % COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
    linrsp(istim, :)  = convCut(stim(istim, :), irf, t_lth);
    numrsp(istim, :)  = linrsp(istim, :).^x.n;
    
    % COMPUTE THE NORMALIZATION DENOMINATOR -------------------------------
    poolrsp(istim, :) = convCut(linrsp(istim, :), irf_norm, t_lth);
    demrsp(istim, :)  = x.sigma.^x.n + poolrsp(istim, :).^x.n;
    
    % COMPUTE HTE NORMALIZATION RESPONSE
    normrsp(istim, :) = x.scale.*(numrsp(istim, :)./demrsp(istim, :));
    
%     normrsp(istim, :) = normRange(normrsp(istim, :));
     
    if ismember(istim, round((1:10)/10* size(stim, 1)-1)) % every 10% draw a dot
        fprintf('[%s]: (%2.0f%%/irf)\n',mfilename,100*(istim/size(stim, 1)));drawnow;
    end
end
result{1} = normrsp';
% result{2} = linrsp';
% result{3} = numrsp';
% result{4} = demrsp';


end