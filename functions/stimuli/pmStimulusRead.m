function stimValues = pmStimulusRead(input)
% Check if it is a string, if it is reads the values and return values
%
%  Input:
%    input:  Either a path to a file or the stimulus movie
%    (row,col,time).
%
%  Output: values
%    stimValues:  The output values
%
% See also
%

if (isstring(input) || ischar(input))
    if exist(input,'file')
        % File must contain the variable 'stim' which includes values
        % and ....
        s = load(input);
        stimValues = s.stim;
    else
        error('File not found %s\n',input);
    end
else
    stimValues = input;
end


end

