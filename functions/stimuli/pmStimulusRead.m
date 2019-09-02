function stimValues = pmStimulusRead(input, varargin)
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


% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired ('input');
p.addParameter('format','matlab', @ischar);
p.parse(input,varargin{:});
format = p.Results.format;


if (isstring(input) || ischar(input))  
    % TODO: make the system write relative paths... 
    % Rigth now in the datatables travel absolute paths, and this has
    % problems when reading local data (stim files as here).
    % For now do a check and it the root it not the same as pmRootPath,
    % change it.
    [rPath,fn,ext] = fileparts(input);
    rootPathParts  = split(rPath,filesep);
    prfmodelLoc    = find(ismember(rootPathParts,'PRFmodel'),1);
    rootPath       = strjoin(rootPathParts(1:prfmodelLoc)',filesep);
    restPath       = strjoin(rootPathParts((prfmodelLoc+1):end)',filesep);
    
    if ~(strcmp(rootPath,pmRootPath))
        input = fullfile(pmRootPath,restPath,[fn ext]);
    end
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

stimValues = double(stimValues);

switch format
    case 'matlab'
        
    case 'nifti'
        [r,c,t] = size(stimValues);
        stimValues = reshape(stimValues, [r,c,1,t]);
        
    otherwise
        error('Unkown format %s',format)
end



end

