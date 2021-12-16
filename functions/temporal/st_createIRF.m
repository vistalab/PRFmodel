function tmodel= st_createIRF(pm)
% pm=prfModel;
if (~exist('synBOLD', 'dir')); mkdir('synBOLD'); end%if
irf_file = [pm.Temporal.IRFpath pm.Temporal.Name '_irf.mat'];

% load default and preprocess accordingly temporal design
temp_type = char(pm.Temporal.temporalModel);
seqType = char(pm.Temporal.stimseq);

stimValues = pm.Stimulus.getStimValues;

[reps, trial_dur,blank_dur]= st_stimrep(stimValues);
temporal.trial_dur = trial_dur;
temporal.blank_dur = blank_dur;

switch seqType 
    case {'a'}
        temporal.stim_on = temporal.trial_dur;
        temporal.stim_off = 0;
        temporal.stim_num = inf;
    case {'b'}
        temporal.stim_on = 1/30;
        temporal.stim_off = 4/30;
        temporal.stim_num = temporal.trial_dur / (temporal.stim_on + temporal.stim_off);
    case {'c'}
        temporal.stim_on = 0.3;
        temporal.stim_off = 1/30;
        temporal.stim_num = temporal.trial_dur / (temporal.stim_on + temporal.stim_off);
       
    otherwise
        temporal.stim_on = str2num(pm.Temporal.stim_on)*1/str2num(pm.Temporal.resolution);
        temporal.stim_off = str2num(pm.Temporal.stim_off)*1/str2num(pm.Temporal.resolution);
        temporal.stim_num = temporal.trial_dur / (temporal.stim_on + temporal.stim_off);
        temporal.stim_num = round(temporal.stim_num, 6); % for decimal errors
end

if nearZero(rem(temporal.stim_num,1)) == 1
    disp(temporal)
elseif temporal.stim_num == inf
    disp(temporal)
else
     error('Stimulus trial is not a whole number')
end

% st = st_spaceTime_IRF(stimValues,temporal);
st = st_stimconvert2(stimValues,temporal);
%     [st] = st_stimconvert(stimValues, ...
%         char(pm.Temporal.stimseq), ...
%         0.033 ,30);

% 
% % convert to ms with different experiment temporal types
% if contains(irf_file,'abc')
%     [st1] = st_stimconvert(stimValues, ...
%         char('a'), ...
%         0.033 ,30);
%     [st2] = st_stimconvert(stimValues, ...
%         char('b'), ...
%         0.033 ,30);
%     [st3] = st_stimconvert(stimValues, ...
%         char('c'), ...
%         0.033 ,30);
%     st.stim = cat(3, st1.stim, st2.stim, st3.stim);
% else
%     [st] = st_stimconvert(stimValues, ...
%         char(pm.Temporal.stimseq), ...
%         0.033 ,30);
% end
 
    
ststim = reshape(st.stim, ...
    size(st.stim,1)*size(st.stim,2),[]);

allstimimages = ststim';
cellimage = num2cell(allstimimages,1)';

mn = length(cellimage);
ms = [[1:ceil(mn./200):mn-2] mn+1];

stimirf_chan_s = zeros(size(allstimimages),'single');
stimirf_chan_t = zeros(size(allstimimages),'single');
stimirf = zeros(size(allstimimages),'single');
fprintf('Generating irf for %s model...\n', temp_type)

for mn = 1:numel(ms)-1

    tmodel = stModel(temp_type,{'Exp2'},'default');
    tmodel.stim = cellimage(ms(mn):ms(mn+1)-1);
    
    [tmodel.onsets, tmodel.offsets, dur] = cellfun(@st_codestim, ...
        tmodel.stim, 'uni', false);
    
    stimirf_base = zeros(size(ststim,2),size(tmodel.stim,1));
    stimirf_t_s = zeros(size(ststim,2),size(tmodel.stim,1));
    stimirf_t_t = zeros(size(ststim,2),size(tmodel.stim,1));
    
    %%%% deal with empty cells (all zero)
    empty_cells = cellfun(@isempty, tmodel.onsets);
    if sum(empty_cells) > 0
        sample =  tmodel.stim(find(empty_cells,1));
        if tmodel.num_channels == 2
            emptysample = {[single(sample{1}) single(sample{1})]};
        else
            emptysample = {[single(sample{1})]};
        end
        tmodel.stim(empty_cells) ={[]};
    end
    
    %%%% make IRFs
    dohrf = 2;
    tmodel = st_pred_runs(tmodel,dohrf);
    
    %%% put dummy into empty cells
    if sum(empty_cells) > 0
        tmodel.pixel_preds(empty_cells) = emptysample;
    end
    
    %%%% transpose and sum up two channels
    stimirf_base = cellfun(@transpose,tmodel.pixel_preds,'UniformOutput',false);
    
    if contains(temp_type,'2ch')
        stimirf_t_s  = cell2mat(cellfun(@(X) X(1,:), stimirf_base, 'uni', false))'; %  sustained
        stimirf_t_t  = cell2mat(cellfun(@(X) X(2,:), stimirf_base, 'uni', false))'; % transient
        stimirf_base = cell2mat(cellfun(@sum, stimirf_base, 'uni', false))';        % sum
        
        stimirf_chan_s(:,ms(mn):ms(mn+1)-1) = stimirf_t_s;
        stimirf_chan_t(:,ms(mn):ms(mn+1)-1) = stimirf_t_t;
    else
        stimirf_t_s  = cell2mat(cellfun(@(X) X(1,:), stimirf_base, 'uni', false))'; %  sustained
        stimirf_chan_s(:,ms(mn):ms(mn+1)-1) = stimirf_t_s;
    end
    
    %                         stimirf(:,ms(mn):ms(mn+1)-1) = stimirf_base;
    
    if ismember(mn, round((1:10)/10* numel(ms)-1)), % every 10% draw a dot
        fprintf(1,'(irf)');drawnow;
    end
    
    
end
% made single to account for matrix storage size
if tmodel.num_channels == 2
    tmodel.chan_preds{1} = stimirf_chan_s;
    tmodel.chan_preds{2} = stimirf_chan_t;
else
    tmodel.chan_preds{1} = stimirf_chan_s;
end


% temporal channel normalization
% need to think about ways of normalizing in the future.
temporal_channel_normalization=false;
if temporal_channel_normalization
    tmodel.normT = max(max(stimirf_chan_s))/max(max(stimirf_chan_t));
else
    tmodel.normT = 1;
end


save(irf_file, 'tmodel','-v7.3');

end
