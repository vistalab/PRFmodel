function model = pred_trials_1ch_glm(model)
% Generates trial predictors using a simple general linear model.

% get design parameters
sessions = model.sessions; nsess = length(sessions); irfs = model.irfs;
cond_list = model.cond_list; nconds_max = max(cellfun(@length, cond_list));
fs = model.fs; tr = model.tr; nexps = model.num_exps;
model.trial_preds.pred = cell(nconds_max, nsess, nexps);
stimfiles = model.stimfiles; nruns = model.num_runs; rcnt = 1;

% generate trial predictors for each session
for ee = 1:nexps
    [~, ~, ~, ~, ton, toff, tc, ~, cl] = tch_stimfile(stimfiles{rcnt, 1});
    for cc = 1:length(cond_list{ee})
        % find trial onset and offset times and calculate duration
        idx = find(strcmp(cl{cc}, tc), 1);
        td = max(tr, ceil(toff(idx) - .001) - ton(idx));
        % extract stimulus vector from condition time window
        cstim_start = round(fs * (ton(idx) - model.pre_dur)) + 1;
        cstim_stop = round(fs * (ton(idx) + td + model.post_dur));
        cstim = model.stim{rcnt, 1}(cstim_start:cstim_stop, :);
        cstim(1:fs * model.pre_dur, :) = 0;
        cstim(fs * (model.pre_dur + td):size(cstim, 1), :) = 0;
        for ss = 1:length(sessions)
            stim_trs = reshape(cstim', size(cstim, 2), tr * fs, []);
            pred = ceil(squeeze(mean(stim_trs, 2))');
            fmri = convolve_vecs(pred, canonical_hrf(tr), 1, 1);
            model.trial_preds.pred{cc, ss, ee} = fmri;
        end
    end
    rcnt = rcnt + nruns(ee, 1);
end

end
