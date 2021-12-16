function all_sessions = find_sessions(project_dir)
% Finds paths to all sessions in project data directory and returns an
% error if none are found.
% 
% INPUT: path to project directory
% OUTPUT: paths to all session directories in project data directory
% 
% AS 5/2017
%
% edited to include full path
% KIS 2/2020


data_dirs = dir(fullfile(project_dir, 'data'));
fulldir = {fullfile(project_dir, 'data/')};

all_sessions = {data_dirs([data_dirs.isdir]).name};
all_sessions(ismember(all_sessions, {'.' '..'})) = [];

fulldir = repmat(fulldir,1,length(all_sessions));

all_sessions = cellfun(@(c,x) strcat(c,x) ,fulldir,all_sessions,'UniformOutput',false);

if isempty(all_sessions)
    error('No sessions found in data directory.');
end
    
end

