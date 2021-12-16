function model = pred_runs_2ch_lin_quad(model,dohrf)
% Generates run predictors using the 2 temporal-channel with optimized 
% linear sustained and quadratic transient channels. 
if nargin< 2
    dohrf = 1;
end

% get design parameters
fs = model.fs; tr = model.tr; stim = model.stim;
nruns_max = size(stim, 1); empty_cells = cellfun(@isempty, stim);
irfs_names = fieldnames(model.irfs); irfs = [];
for ff = 1:length(irfs_names)
    iname = model.irfs.(irfs_names{ff});
    irfs.(irfs_names{ff}) = repmat(iname, nruns_max, 1);
end

% generate run predictors for each session
predS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    stim, irfs.nrfS, 'uni', false); predS(empty_cells) = {[]};
predTq = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs) .^ 2, ...
    stim, irfs.nrfT, 'uni', false); predTq(empty_cells) = {[]};

if dohrf == 1
    fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
        predS, irfs.hrf, 'uni', false); fmriS(empty_cells) = {[]};
    fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
        predTq, irfs.hrf, 'uni', false); fmriT(empty_cells) = {[]};
    run_preds = cellfun(@(X, Y) [X Y * model.normT], fmriS, fmriT, 'uni', false);
    
    model.run_preds = run_preds;
else
    model.normT = 1;
    model.pixel_preds = cellfun(@(X, Y) [X Y * model.normT], ...
        predS, predTq, 'uni', false);
end



end