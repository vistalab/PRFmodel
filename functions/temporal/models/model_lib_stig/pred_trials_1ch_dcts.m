function model = pred_trials_1ch_dcts(model)
% Generates trial predictors using the single-channel dCTS model proposed 
% by Zhou et al. (2017). 

% get design parameters
sessions = model.sessions; nsess = length(sessions); irfs = model.irfs;
cond_list = model.cond_list; nconds_max = max(cellfun(@length, cond_list));
fs = model.fs; tr = model.tr; nexps = model.num_exps;
model.trial_preds.pred = cell(nconds_max, nsess, nexps);
stimfiles = model.stimfiles; nruns = model.num_runs; rcnt = 1;

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
        for ss = 1:length(sessions)
            % convolve stimulus with channel IRF and lpf
            predSl = convolve_vecs(cstim, irfs.nrfS{ss}, fs, fs);
            predSd = convolve_vecs(predSl, irfs.lpf{ss}, fs, fs);
            predS = (predSl .^ 2) ./ (model.params.sigma{ss} .^ 2 + predSd .^ 2);
            % convolve neural predictors with HRF
            fmri = convolve_vecs(predS, irfs.hrf{ss}, fs, 1 / tr);
            model.trial_preds.pred{cc, ss, ee} = fmri;
        end
    end
    rcnt = rcnt + nruns(ee, 1);
end

end
