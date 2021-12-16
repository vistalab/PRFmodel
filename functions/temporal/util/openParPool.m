function openParPool(nCPU)
% first check if sherlock CPU enviornment is there.
% if not the function opens a given number of cpus

% default is 16
if notDefined('nCPU')
    nCPU = 4;
end

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    pc = parcluster('local');
    if ~isempty(getenv('SLURM_CPUS_ON_NODE'))
%         parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
        parpool(pc, nCPU);
    else
        parpool(pc, nCPU);
        fprintf('default number of parpool nodes loaded \n');
    end
end

end