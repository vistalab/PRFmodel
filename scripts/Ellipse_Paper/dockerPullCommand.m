function dockerPullCommand(image)
%DOCKERPULLCOMMAND Checks if the required image is downloaded
    [s,r] = system(sprintf('docker pull %s',image));
    if s > 0
        fprintf('\n %s not installed or Docker client not running. Error is:\n',image)
        fprintf(' --> %s',r)
    else
        fprintf('\n %s installed succesfully\n',image)
    end
end

