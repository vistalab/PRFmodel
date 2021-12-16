function st_saveBOLD(DT,output_dir)

%%
% if (~exist('synBOLD', 'dir')); mkdir('synBOLD'); end%if

% fname = DT.pm(1).Temporal.Name; 
% RFs = [DT.pm(1).RF.Centerx0 DT.pm(1).RF.Centery0  DT.pm(1).RF.sigmaMajor];
% RFs = num2str(RFs);
% RFs(RFs == ' ') = [];
% 
fname1 = [output_dir '/synBOLD.mat'];
fname2 = [output_dir '/table.mat'];
for ee = 1:height(DT)
    chan{ee}    = cell2mat(DT(ee,:).pm.Temporal.chan_preds)';
    synBOLD{ee} = DT(ee,:).pm.BOLDnoise';
end

save(fname1,'chan','synBOLD')
save(fname2,'DT')



end