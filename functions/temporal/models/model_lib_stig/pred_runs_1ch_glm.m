% 3/2020 edited by KIS
% fixed the dimension problem in line 21

function model = pred_runs_1ch_glm(model,dohrf)
% Generates run predictors using a simple general linear model.
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
stim_trs = cellfun(@(X) reshape(X', size(X, 2), tr * fs, []), ...
    stim, 'uni', false); stim_trs(empty_cells) = {[]};

stim_trs = cellfun(@(X) ceil(squeeze(mean(X, 2))'), ...
    stim_trs, 'uni', false); stim_trs(empty_cells) = {[]};



if dohrf == 1
    stim_trs = cellfun(@(X) reshape(X', size(X, 2), tr * fs, []), ...
        stim, 'uni', false); stim_trs(empty_cells) = {[]};
    
    stim_trs = cellfun(@(X) ceil(squeeze(mean(X, 2))'), ...
        stim_trs, 'uni', false); stim_trs(empty_cells) = {[]};
    
    % fixed bug to work both for PNAS and Plosone dataset
    if size(stim_trs{1},1) <  size(stim_trs{1},2)
        stim_trs = cellfun(@transpose,stim_trs,'UniformOutput',false);
    end
    
    run_preds = cellfun(@(X) convolve_vecs(X, canonical_hrf(tr), 1, 1), ...
        stim_trs, 'uni', false); run_preds(empty_cells) = {[]};
    model.run_preds = run_preds;
    
    
else
    
    model.pixel_preds = stim;
    
end


end


