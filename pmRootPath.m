function rootPath = pmRootPath()
% Determine path to root of the PRFmodel directory
%
%        rootPath = pmRootPath;
%
% This function MUST reside in the directory at the base of the
% PRFmodel directory structure 
%
% Copyright Stanford team, mrVista, 2019
% test
rootPath = which('pmRootPath');

rootPath = fileparts(rootPath);

rootPath = '/export/home/glerma/testFolder/PRFmodel/ellipse'

return
