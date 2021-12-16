function model = pred_runs_2ch_lin_htd(model)
% Generates run predictors using the 2-channel hemodynamic temporal 
% derivative model proposed by Henson et al. (2002). 

% get design parameters
fs = model.fs; tr = model.tr; stim = model.stim;
nruns_max = size(stim, 1); empty_cells = cellfun(@isempty, stim);
irfs_names = fieldnames(model.irfs); irfs = [];
for ff = 1:length(irfs_names)
    iname = model.irfs.(irfs_names{ff});
    irfs.(irfs_names{ff}) = repmat(iname, nruns_max, 1);
end

% generate run predictors for each session
fmri_c = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    stim, irfs.hrf, 'uni', false); fmri_c(empty_cells) = {[]};
fmri_d = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    stim, irfs.dhrf, 'uni', false); fmri_d(empty_cells) = {[]};
run_preds = cellfun(@(X, Y) [X Y], fmri_c, fmri_d, 'uni', false);
model.run_preds = run_preds;

end
