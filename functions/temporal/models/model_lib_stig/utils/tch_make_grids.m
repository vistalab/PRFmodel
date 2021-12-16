function tch_make_grids(model_type, exps)

roi = tchROI('pFus_faces', exps);
roi = select_sessions(roi);
sessions = roi.sessions;

for ss = 1:length(sessions)
    roi = tchROI('pFus_faces', exps, sessions{ss});
    roi = tch_runs(roi);
    model = tchModel(model_type, roi.experiments, sessions{ss});
    model = code_stim(model);
    model = pred_runs(model);
    model = pred_trials(model);
    roi = tch_trials(roi, model);
    [roi, model] = tch_fit(roi, model);
    [rois, models] = tch_grid_search(roi, model, 1, 5);
end

end
