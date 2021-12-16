% st_tModel: Code for modeling fMRI responses to time-varying visual stimuli.
%
% CONSTRUCTOR INPUTS
%   1) type: which model to use
%      Hemodynamic models:
%        '1ch-glm'     -- general linear model for fMRI (Boynton 1996)
%        '2ch-lin-htd' -- hemodynamic temporal derivative (HTD; Henson 2002)
%        '1ch-balloon' -- nonlinear balloon model (Buxton, 1998)
%      Optimized single-channel models:
%        '1ch-lin'  -- linear sustained channel with neural IRF
%        '1ch-exp'  -- susatined channel with adaptation (A)
%        '1ch-pow'  -- CTS with power law (CTS-p; Zhou 2017)
%        '1ch-div'  -- CTS with divisive normalization (CTS-n; Zhou 2017)
%        '1ch-dcts' -- dynamic CTS (dCTS; Zhou 2017)
%        '1ch-rect' -- transient channel with rectification
%        '1ch-quad' -- transient channel with quadratic nonlinearity (X^2)
%        '1ch-sig'  -- transient channel with sigmoid nonlinearities
%      Optimized two-channel models:
%        '2ch-lin-quad'  -- linear susatined and quadratic transient (L+Q)
%        '2ch-lin-rect'  -- linear susatined and rectified transient (L+R)
%        '2ch-exp-quad'  -- adapted sustained and quadratic transient (A+Q)
%        '2ch-exp-rect'  -- adapted sustained and rectified transient (A+R)
%        '2ch-pow-quad'  -- sustained with CTS-p and quadratic transient (C+Q)
%        '2ch-pow-rect'  -- sustained with CTS-p and rectified transient (C+R)
%        '2ch-lin-sig'   -- linear sustained and sigmoid transient (L+S)
%        '2ch-exp-sig'   -- adapted sustained and sigmoid transient (A+S)
%

%   roi = select_sessions(roi);
%   model = stModel('2ch-lin-quad', exps, roi.sessions);


function result = st_tModel(tModel, param, stim, time)

% fprintf('Generating irf for %s model...\n', tModel)

switch tModel
    case {'glm','1ch-glm'}
        result  = linearModel(param, stim, time);
    case {'DN','1ch-dcts'}
        result = DNmodel(param, stim, time);
    case {'2ch','2ch-exp-sig'}
        result  = twoChansmodel(param, stim, time);
    case {'2ch-lin-sig'}
        result  = twoChansmodel_lin_sig(param, stim, time);
    case {'2ch-css-sig'}
        result  = twoChansmodel_Scss_Tsig(param, stim, time);
        
        %     case '1ch-glm'
        %         model = pred_runs_1ch_glm(model);
        %     case '2ch-lin-htd'
        %         model = pred_runs_2ch_lin_htd(model);
        %     case '1ch-balloon'
        %         model = pred_runs_1ch_balloon(model);
        %     case '1ch-lin'
        %         model = pred_runs_1ch_lin(model);
        %     case '1ch-exp'
        %         model = pred_runs_1ch_exp(model);
        %     case '1ch-pow'
        %         model = pred_runs_1ch_pow(model);
        %     case '1ch-div'
        %         model = pred_runs_1ch_div(model);
        %     case '1ch-rect'
        %         model = pred_runs_1ch_rect(model);
        %     case '1ch-quad'
        %         model = pred_runs_1ch_quad(model);
        %     case '1ch-sig'
        %         model = pred_runs_1ch_sig(model);
        %     case '2ch-lin-rect'
        %         model = pred_runs_2ch_lin_rect(model);
        %     case '2ch-exp-rect'
        %         model = pred_runs_2ch_exp_rect(model);
        %     case '2ch-lin-quad'
        %         model = pred_runs_2ch_lin_quad(model);
        %     case '2ch-exp-quad'
        %         model = pred_runs_2ch_exp_quad(model);
        %     case '2ch-pow-quad'
        %         model = pred_runs_2ch_pow_quad(model);
        %     case '2ch-pow-rect'
        %         model = pred_runs_2ch_pow_rect(model);
        %     case '2ch-lin-sig'
        %         model = pred_runs_2ch_lin_sig(model);
end



end
