function [data_tmp_roi,vert]=pmMaskFsaverage(data_tmp,projectDir,subject,roiname)

data_tmp_roi = [];

fspth = fullfile(projectDir, 'derivatives', 'freesurfer', ['sub-' subject]);
%     path2roi = {'V1_exvivo';'V2_exvivo'};


hemispheres = {'lh';'rh'};
lcurv  = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
rcurv  = read_curv(fullfile(fspth, 'surf', 'rh.curv'));
idx{1} = 1:numel(lcurv);
idx{2} = (1:numel(rcurv))+numel(lcurv);

for hemi = 1 : length(hemispheres)
    
    
    
    roi = [];
    
    for r = 1 : length(roiname)
        
        ind  = read_label(['sub-' subject],sprintf ('%s.%s%s',hemispheres{hemi},roiname{r}));
        roi  = [roi; ind(:,1) + 1];
        
    end
    
    data_hemi = data_tmp(idx{hemi},:,:,:);
    data_hemi_roi = data_hemi(roi,:,:,:);
    
    data_tmp_roi = cat(1,data_tmp_roi,data_hemi_roi);
    vert{hemi,:}  = roi;
    
    
end

end

