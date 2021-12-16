function model = pred_runs_2ch_pow_quad(model)
% Generates run predictors using the 2 temporal-channel model with CTS-pow
% on sustained and quadratic transient channel. 

% get design parameters
fs = model.fs; tr = model.tr; stim = model.stim;
nruns_max = size(stim, 1); empty_cells = cellfun(@isempty, stim);
params_names = fieldnames(model.params); params = [];
for pp = 1:length(params_names)
    pname = model.params.(params_names{pp});
    params.(params_names{pp}) = repmat(pname, nruns_max, 1);
end
irfs_names = fieldnames(model.irfs); irfs = [];
for ff = 1:length(irfs_names)
    iname = model.irfs.(irfs_names{ff});
    irfs.(irfs_names{ff}) = repmat(iname, nruns_max, 1);
end

% generate run predictors for each session
predSp = cellfun(@(X, Y, Z) convolve_vecs(X, Y, fs, fs) .^ Z, ...
    stim, irfs.nrfS, params.epsilon, 'uni', false); predSp(empty_cells) = {[]};
predTq = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs) .^ 2, ...
    stim, irfs.nrfT, 'uni', false); predTq(empty_cells) = {[]};
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predSp, irfs.hrf, 'uni', false); fmriS(empty_cells) = {[]};
fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predTq, irfs.hrf, 'uni', false); fmriT(empty_cells) = {[]};
run_preds = cellfun(@(X, Y) [X Y * model.normT], fmriS, fmriT, 'uni', false);
model.run_preds = run_preds;

end
