function [status,result] = dockerPullCommand(image)
%DOCKERPULLCOMMAND Checks if the required image is downloaded
%
% Synopsis
%   dockerPullCommand(image)
%
% Input:
%   image -   The path and name of the Docker image  (string).  
%      For example,  'garikoitz/prfsynth:2.0.0'
%      
% Output:
%    status -
%    result -
%
% See also

[status,result] = system(sprintf('docker pull %s',image));

if status > 0
    fprintf('\n %s not installed or Docker client not running. Error is:\n',image)
    fprintf(' --> %s',result)
else
    fprintf('\n %s installed succesfully\n',image)
end

end

