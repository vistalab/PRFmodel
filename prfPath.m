function rootPath = prfPath()
% Determine path to root of the mrVista directory
%
%        rootPath = vistaRootPath;
%
% This function MUST reside in the directory at the base of the
% afqDimensionality directory structure 
%
% Copyright Stanford team, mrVista, 2018

rootPath = which('prfPath');

rootPath = fileparts(rootPath);

return
