function yesno = isclose(x,y,varargin)
%isclose: when isequal gives 0 but you know the vectors are close within a
%         tolerance
%   x,y       = numbers or vectors
%   tolerance = the difference between numbers should be less than


%% Read the inputs

% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired ('x'            ,         @isnumeric);
p.addRequired ('y'            ,         @isnumeric);
p.addParameter('tolerance'    , 10^-10, @isnumeric);
p.addParameter('returnvector' , false , @islogical);

p.parse(x, y, varargin{:});

tolerance       = p.Results.tolerance;
returnvector    = p.Results.returnvector;

%% Do the thing
    if returnvector
        if isequal(size(x),size(y))
            yesno = abs(x-y) < tolerance;
        else
            yesno = false;
        end
    else
        yesno = isequal(size(x),size(y)) ...
            && ...
            all(abs(x(:)-y(:)) < tolerance);
    end
end






