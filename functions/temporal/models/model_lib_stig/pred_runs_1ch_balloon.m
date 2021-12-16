function model = pred_runs_1ch_balloon(model)
% Generates run predictors using a standard single-channel balloon model.

% define global variables
global which_tau tauN tauP tauMTT alpha E0 V0 tau;
% get design parameters
params = model.params; irfs_init = model.irfs; dt = params.delta_t;
fs = model.fs; tr = model.tr; rd = model.run_durs; stim = model.stim;
cat_list = unique([model.cats{:}]); ncats = length(cat_list);
[nruns_max, nsess] = size(model.run_durs); empty_cells = cellfun(@isempty, rd);
irfs_names = fieldnames(model.irfs); irfs = [];
for ff = 1:length(irfs_names)
    iname = model.irfs.(irfs_names{ff});
    irfs.(irfs_names{ff}) = repmat(iname, nruns_max, 1);
end
time_vecs = cellfun(@(X) [0:dt:X]', rd, 'uni', false);

%% model nonlinear hemodynamics with balloon model

% convolve input with a gamma and normalize
in_flow = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    stim, irfs.gamma, 'uni', false);
in_flow = cellfun(@(X, Y) 0.7 * (X / sum(Y)) + 1, ...
    in_flow, irfs.gamma, 'uni', false);

% get the simulated values of all variables for each predictor in each run
run_preds = cell(nruns_max, nsess);
for ss = 1:nsess
    [~, session_id] = fileparts(model.sessions{ss});
    fprintf('Simulating balloon model for %s', session_id);
    fname = ['1ch_balloon_model_' [model.experiments{:}] '.mat'];
    fpath = fullfile(model.sessions{ss}, 'Stimuli', fname);
    if exist(fpath, 'file') == 2
        load(fpath);
        run_preds(:, ss) = session_run_preds;
        fprintf('\n');
    else
        for rr = 1:nruns_max
            fprintf('.');
            for pp = 1:size(in_flow{rr, ss}, 2)
                % initialize variables
                [v, q, IN_FLOW, OUT_FLOW, CMRO2] = deal(1);
                S = 0; OEF = params.E0; t = time_vecs{rr, ss}; 
                which_tau = 1; tauMTT = params.tauMTT; tau = params.tauP;
                tauP = params.tauP; tauN = params.tauN;
                alpha = params.alpha; E0 = params.E0; V0 = params.V0;
                % get the simulated values of all variables
                for ii = 1:size(time_vecs{rr, ss}, 1) - 1
                    ii_flow = in_flow{rr, ss}(ii, pp);
                    v(ii + 1) = runge_kutta(dt, @dvdt, t(ii), v(ii), ii_flow);
                    OUT_FLOW(ii + 1) = flow_out(v(ii + 1), t(ii), ii_flow);
                    q(ii + 1) = runge_kutta(dt, @dqdt, t(ii), q(ii), v(ii), ii_flow);  
                    S1 = params.k1 * (1 - q(ii + 1));
                    S2 = params.k2 * (1 - q(ii + 1) / v(ii + 1));
                    S3 = params.k3 * (1 - v(ii + 1));
                    S(ii + 1) = V0 * (S1 + S2 + S3);
                    OEF(ii + 1) = 1 - (1 - E0) .^ (1 ./ ii_flow);
                    CMRO2(ii + 1) = (OEF(ii + 1) / E0) * IN_FLOW(ii);
                    IN_FLOW(ii + 1) = ii_flow;
                end
                Sr = convolve_vecs(S(1:rd{rr, ss} / params.delta_t)', 1, fs, 1 / tr);
                run_preds{rr, ss}(:, pp) = Sr;
            end
        end
        session_run_preds = run_preds(:, ss); 
        save(fpath, 'session_run_preds', '-v7.3');
        % store predictors for each experiment separately as welll
        if length(model.experiments) > 1
            srp = session_run_preds; rcnt = 0;
            for ee = 1:length(model.experiments)
                session_run_preds = srp(rcnt + 1:rcnt + model.num_runs(ee, ss));
                fname = ['balloon_model_' model.experiments{ee} '.mat'];
                fpath = fullfile(model.sessions{ss}, 'Stimuli', fname);
                save(fpath, 'session_run_preds', '-v7.3');
                rcnt = rcnt + model.num_runs(ee, ss);
            end
        end
        fprintf('\n');
    end
end
run_preds(empty_cells) = {[]};
model.run_preds = run_preds;

end
