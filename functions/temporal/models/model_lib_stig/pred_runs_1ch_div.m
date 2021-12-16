function model = pred_runs_1ch_div(model,dohrf)
% Generates run predictors using the single-channel CTS-div model proposed 
% by Zhou et al. (2017).

if nargin< 2
    dohrf = 1;
end

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
predSl = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    stim, irfs.nrfS, 'uni', false);
predSn = cellfun(@(X) X .^ 2, predSl, 'uni', false);
predSd = cellfun(@(X, Y) X + Y .^ 2, predSn, params.sigma, 'uni', false);
predS = cellfun(@(X, Y) X ./ Y, predSn, predSd, 'uni', false);

if dohrf == 1
    run_preds = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
        predS, irfs.hrf, 'uni', false); run_preds(empty_cells) = {[]};
    model.run_preds = run_preds;
else
    model.pixel_preds = predS;
    model.pixel_preds = cellfun(@single,model.pixel_preds,'un',0);
    
end




end
