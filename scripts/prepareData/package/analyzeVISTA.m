function analyzeVISTA(mainDir,cfg,averageFolDir,subject,session,apertureFolder,dockerscript,scriptDir)



stimulusDir = [apertureFolder filesep sprintf('sub-%s/ses-%s/',subject,session)];
cfgfile = [cfg.average_filename '_cfg.json'];

system(sprintf('chmod u+x %s/%s',scriptDir,dockerscript))


% if contains(dockerscript,'singularity')
    
tmp = [tempdir 'config_pRF'];
mkdir(tmp)
% copyfile(cfgfile,[tmp filesep cfg.param.basename '_cfg.json'])
copyfile(cfgfile,[tmp filesep 'config.json'])



if contains(dockerscript,'singularity')
    
    [~,user] = system('echo $USER');
    singimg = sprintf('/home/%s/vista.simg',user(1:end-1));
    system(sprintf('%s/%s vista %s %s %s %s %s',scriptDir,dockerscript,averageFolDir,tmp,stimulusDir,mainDir,singimg));
    
elseif contains(dockerscript,'docker')
    
    system(sprintf('%s/%s vista %s %s %s %s',scriptDir,dockerscript,averageFolDir,tmp,stimulusDir,mainDir));
    
end

