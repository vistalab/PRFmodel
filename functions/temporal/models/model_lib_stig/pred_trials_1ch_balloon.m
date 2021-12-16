function model = pred_trials_1ch_balloon(model)
% Generates trial predictors using a standard single-channel balloon model.

% define global variables and signal parameters
global which_tau tauN tauP tau tauMTT alpha E0 V0;

% get design parameters
sessions = model.sessions; nsess = length(sessions); irfs = model.irfs;
cond_list = model.cond_list; nconds_max = max(cellfun(@length, cond_list));
fs = model.fs; tr = model.tr; nexps = model.num_exps; dt = model.params.delta_t;
model.trial_preds.pred = cell(nconds_max, nsess, nexps);
stimfiles = model.stimfiles; nruns = model.num_runs; rcnt = 1;

% generate trial predictors for each session
for ee = 1:nexps
    [~, ~, ~, ~, ton, toff, tc, ~, cl] = tch_stimfile(stimfiles{rcnt, 1});
    for cc = 1:length(cond_list{ee})
        % find trial onset and offset times and calculate duration
        idx = find(strcmp(cl{cc}, tc), 1);
        td = ceil(toff(idx) - .001) - ton(idx);
        % extract stimulus vector from condition time window
        cstim_start = round(fs * (ton(idx) - model.pre_dur)) + 1;
        cstim_stop = round(fs * (ton(idx) + td + model.post_dur));
        cstim = model.stim{rcnt, 1}(cstim_start:cstim_stop, :);
        cstim(1:fs * model.pre_dur, :) = 0;
        cstim(fs * (model.pre_dur + td):size(cstim, 1), :) = 0;
        t = 0:model.params.delta_t:model.pre_dur + td + model.post_dur;
        for ss = 1%:length(sessions)
            in_flow = convolve_vecs(cstim, irfs.gamma{ss}, fs, fs);
            in_flow = 0.7 * (in_flow / sum(irfs.gamma{ss})) + 1;
            for pp = 1:size(cstim, 2)
                % initialize variables
                [v, q, IN_FLOW, OUT_FLOW, CMRO2] = deal(1);
                S = 0; OEF = model.params.E0; which_tau = 1;
                tau = model.params.tauP; tauMTT = model.params.tauMTT;
                tauP = model.params.tauP; tauN = model.params.tauN;
                E0 = model.params.E0; V0 = model.params.V0;
                alpha = model.params.alpha;
                % get the simulated values of all variables                
                for ii = 1:length(t) - 1
                    ii_flow = in_flow(ii, pp);
                    v(ii + 1) = runge_kutta(dt, @dvdt, t(ii), v(ii), ii_flow);
                    OUT_FLOW(ii + 1) = flow_out(v(ii + 1), t(ii), ii_flow);
                    q(ii + 1) = runge_kutta(dt, @dqdt, t(ii), q(ii), v(ii), ii_flow);  
                    S1 = model.params.k1 * (1 - q(ii + 1));
                    S2 = model.params.k2 * (1 - q(ii + 1) / v(ii + 1));
                    S3 = model.params.k3 * (1 - v(ii + 1));
                    S(ii + 1) = V0 * (S1 + S2 + S3);
                    OEF(ii + 1) = 1 - (1 - E0) .^ (1 ./ ii_flow);
                    CMRO2(ii + 1) = (OEF(ii + 1) / E0) * IN_FLOW(ii);
                    IN_FLOW(ii + 1) = ii_flow;
                end
                Sr = convolve_vecs(S(1:length(t) - 1)', 1, fs, 1 / tr);
                model.trial_preds.pred{cc, ss, ee} = Sr;
            end
        end
        for ss = 2:nsess
            model.trial_preds.pred{cc, ss, ee} = model.trial_preds.pred{cc, 1, ee};
        end
    end
    rcnt = rcnt + nruns(ee, 1);
end

end
