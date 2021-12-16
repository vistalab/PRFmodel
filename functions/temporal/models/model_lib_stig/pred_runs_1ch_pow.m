function model = pred_runs_1ch_pow(model)
% Generates run predictors using the single-channel CTS-pow model proposed 
% by Zhou et al. (2017).

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
predS = cellfun(@(X, Y, Z) convolve_vecs(X, Y, fs, fs) .^ Z, ...
    stim, irfs.nrfS, params.epsilon, 'uni', false); predS(empty_cells) = {[]};
run_preds = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predS, irfs.hrf, 'uni', false); run_preds(empty_cells) = {[]};
model.run_preds = run_preds;

end
