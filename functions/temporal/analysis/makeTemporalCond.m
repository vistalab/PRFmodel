function [sequence, trial,ons,offs] = makeTemporalCond()
 %% this works in "st_makePred" file
 load('ExpDesign');

% Run X Time
% amplitude is the temporal cond is the stimulus condition 
% to see result
% imagesc(sequnece)

%% spatial seqeunce
 expseq = [0,0,1:9,0,0,10:18,0,0,19:27,28:36,0,0];
 
 %% get image-wise Seqeunce in ms resolution
ons = [8 8  8  2 12 52 16 48 300]';
offs= [2 12 52 8 8  8  4  12   0]';
ifi = 60;
trial_dur = 5;
stimprm = zeros(length(ons),4);
for i = 1:length(ons)
    stimprm(i,1)=ons(i);
    stimprm(i,2)=offs(i);
    stimprm(i,3)=trial_dur;
    stimprm(i,4)=ifi;
end
[trial,nstim]=st_mkTrial(stimprm,1000);
%multiple to give image number index
for tt = 1:length(trial)
    trial{tt}=trial{tt};
end
 
 
 %% make stimulus seqeunce with blanks & save by each run
blank_trial_dur =5;
temporalMask = expseq~=0;
tmp_seq = zeros(size(temporalMask));
% repelem(expseq,trial_dur);
for run = 1:size(ExpDesign,1)
    tmp = ExpDesign(run,:);
    tmp_seq(temporalMask) = tmp;
    runSeqms = [];
    for cond = 1:length(tmp_seq)
        if tmp_seq(cond) == 0 % blanks
            runSeqms= [runSeqms  repelem(0,blank_trial_dur*1000)];
        elseif tmp_seq(cond) ~= 0  % ons
            tt = trial{tmp_seq(cond)}>0;
            tt = tmp_seq(cond)*tt;
            runSeqms= [runSeqms tt];
        end
    end
    sequence(run,:) = runSeqms;
end

 