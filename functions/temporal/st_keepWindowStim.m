function [keep,msStim] = st_keepWindowStim(pm,keepAllPoints)

if notDefined('keepAllPoints')
    keepAllPoints = 1;
end
% irf_file = [pm.Temporal.IRFpath '/' pm.Temporal.Name '_irf.mat'];


try
    %     tmpSeq = load(fullfile(pm.Stimulus.DataPath, [pm.Stimulus.Name '.mat']));
    tmpSeq = load(pm.Stimulus.values,'sequence','resampled');
    seq = tmpSeq.sequence;
    seq = downsample(seq,1000/pm.Temporal.fs);
    tmt= tmpSeq.resampled;
    
    %% convert to ms Stim
    offMask = zeros([size(tmt,1) size(tmt,2)]);
    msStim = zeros(size(tmt,1),size(tmt,2),length(seq));
    for eachimage = 1:length(seq)
        if seq(eachimage) == 0
            msStim(:,:,eachimage) = offMask;
        elseif seq(eachimage) ~= 0
            msStim(:,:,eachimage) = tmt(:,:,seq(eachimage));
        end
    end
    msStim = reshape(msStim, size(msStim,1)*size(msStim,2),[]);
    
    
    stimwindow = nansum(msStim,2);
    % needs to be changed later on
    if keepAllPoints
        % mark all pixels as pixels to keep
        stimwindow(:) = 1;
        instimwindow = find(stimwindow==1);
    else
        stimwindow = stimwindow > 0;
        stimwindow = stimwindow(:);
        instimwindow = find(stimwindow==1);
    end
    
    keep = instimwindow;
    msStim   = sparse(msStim(keep,:));
catch
%     msStim = pm.Stimulus.getStimValues;
%     msStim = reshape(msStim, size(msStim,1)*size(msStim,2),[]);
%     msStim = resample(msStim,pm.Temporal.fs,1,'Dimension',2);
% 
%     stimwindow = nansum(msStim,2);
%     % needs to be changed later on
%     if keepAllPoints
%         % mark all pixels as pixels to keep
%         stimwindow(:) = 1;
%         instimwindow = find(stimwindow==1);
%     else
%         stimwindow = stimwindow > 0;
%         stimwindow = stimwindow(:);
%         instimwindow = find(stimwindow==1);
%     end
%     
%     keep = instimwindow;
%     msStim   = sparse(msStim(keep,:));
    
    error(1,'no seqeunce struct in the stim file');
end



end