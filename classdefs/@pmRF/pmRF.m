classdef pmRF <   matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for RF-s. Every particular instance of this class
    % will have different parameters, so it will be a children class.
    %
    % Syntax:
    %      rf = RF();
    %
    % Inputs:
    %
    % Outputs:
    %
    % Optional key/value pairs:
    %
    % Description
    %
    % See also
    %
    
    % Examples
    %{
       
    %}
    
    properties
        PM;         % prfModel that has some of the variables we need, such as TR
        Centerx0;   % Deg
        Centery0;   % Deg
        Theta;      % Radians
        sigmaMajor; % Deg
        sigmaMinor; % Deg
        dog_sigmaMajor;
        dog_sigmaMinor;
        dog_Theta;
        dog_Scale; % Scalar, scale factor [0,1]
        Type      ;
        % values;
    end
    properties (SetAccess = private, GetAccess = public)
         values;    % Result. Only can be changes from within this func.
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods (Static)
        function d = defaultsGet
            d.Centerx0       = 0;        % Deg
            d.Centery0       = 0;        % Deg
            d.Theta          = 0;        % Radians
            d.sigmaMajor     = 1;        % Deg
            d.sigmaMinor     = 1;        % Deg
            d.dog_sigmaMajor = 2;
            d.dog_sigmaMinor = 2;
            d.dog_Theta      = 0;
            d.dog_Scale      = 0.5;
            d.Type           = 'mrvista';
            % Convert to table and return
            d = struct2table(d,'AsArray',true);
        end
    end
    methods
        % Constructor
        function rf = pmRF(pm, varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = rf.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('pm'            ,                   @(x)(isa(x,'prfModel')));
            p.addParameter('centerx0'      , d.Centerx0      , @isnumeric);
            p.addParameter('centery0'      , d.Centery0      , @isnumeric);
            p.addParameter('theta'         , d.Theta         , @isnumeric);
            p.addParameter('sigmamajor'    , d.sigmaMajor    , @isnumeric);
            p.addParameter('sigmaminor'    , d.sigmaMinor    , @isnumeric);
            p.addParameter('dog_sigmamajor', d.dog_sigmaMajor, @isnumeric);
            p.addParameter('dog_sigmaminor', d.dog_sigmaMinor, @isnumeric);
            p.addParameter('dog_theta'     , d.dog_Theta     , @isnumeric);
            p.addParameter('dog_scale'     , d.dog_Scale     , @isnumeric);
            p.addParameter('type'          , d.Type{:}       , @ischar);
            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            rf.PM             = pm;
            % The receptive field parameters
            rf.Centerx0       = p.Results.centerx0;
            rf.Centery0       = p.Results.centery0;
            rf.Theta          = p.Results.theta;
            rf.sigmaMajor     = p.Results.sigmamajor;
            rf.sigmaMinor     = p.Results.sigmaminor;
            rf.Type           = p.Results.type;
            rf.dog_sigmaMajor = p.Results.dog_sigmamajor;
            rf.dog_sigmaMinor = p.Results.dog_sigmaminor;
            rf.dog_Theta      = p.Results.dog_theta;
            rf.dog_Scale      = p.Results.dog_scale;
        end
        
        function v = get.TR(rf)
            v      = rf.PM.TR;
        end
        
        % Methods available to this class and childrens, if any
        function compute(rf)
            % Check the sigma values, if we want CSS, we need to multiply the
            % sigma values by the sqrt of the css exponent
            switch rf.PM.Type
                case 'CSS'
                    rf.sigmaMajor     = sqrt(rf.PM.cssexp) * rf.sigmaMajor;
                    rf.sigmaMinor     = sqrt(rf.PM.cssexp) * rf.sigmaMinor;
                    rf.dog_sigmaMajor = sqrt(rf.PM.cssexp) * rf.dog_sigmaMajor;
                    rf.dog_sigmaMinor = sqrt(rf.PM.cssexp) * rf.dog_sigmaMinor;
            end
            % Compute stimulus just in case
            rf.PM.Stimulus.compute;
            % Obtain grid XY
            XY = rf.PM.Stimulus.XY;
            % Calculate values
            if iscell(rf.Type);type = rf.Type{:};else;type=rf.Type;end
            switch type
                case {'mrvista'}
                    rf.values = pmGaussian2d(XY{1}, XY{2}, ...
                                rf.sigmaMajor,rf.sigmaMinor,rf.Theta, ...
                                rf.Centerx0,rf.Centery0);
                case {'analyzeprf'}
                    % res is always square, same side sizes
                    % r and c are not in degrees, are in row and columns, the
                    % code below needs to be modified
                    res     = max(size(XY{1},1), size(XY{1},1));
                    r       = rf.Centery0;
                    c       = rf.Centerx0;
                    sr      = rf.sigmaMajor;
                    sc      = rf.sigmaMinor;
                    xx      = XY{1};
                    yy      = XY{2};
                    ang     = rf.Theta;
                    omitexp = 0;
                    % calculate
                    rf.values = makegaussian2d(res,r,c,sr,sc,xx,yy,ang,omitexp);
                case {'dog'}
                    X        = XY{1};
                    Y        = XY{2};
                    x0       = rf.Centerx0;
                    y0       = rf.Centery0;
                    smaj     = rf.sigmaMajor;
                    smin     = rf.sigmaMinor;
                    theta    = rf.Theta;
                    dogsmaj  = rf.dog_sigmaMajor;
                    dogsmin  = rf.dog_sigmaMinor;
                    dogtheta = rf.dog_Theta;
                    scale    = rf.dog_Scale;
                    % Checks variables
                    % For now we only allow circular DoG, change sigmas and theta = 0
                    smin     = smaj;
                    dogsmin  = dogsmaj;
                    theta    = 0;
                    dogtheta = 0;
                    % Assert smaj < dogsmaj
                    assert(smaj < dogsmaj,'sigma of the internal RF has to be smaller than the external one')
                    % Scale value needs to be [0,1]
                    assert((scale>0) && (scale<1),'DoG scale is %f (should be [0,1])',scale)
                    
                    % Create the gaussians and normalize them (we want same areas in 2D gaussians)
                    rf1 = pmGaussian2d(X,Y,smaj,smin,theta,x0,y0);
                    rf1 = rf1 ./ sqrt(smaj.*2.*pi.*smin);
                    rf2 = pmGaussian2d(X,Y,dogsmaj,dogsmin,dogtheta,x0,y0);
                    rf2 = rf2 ./ sqrt(dogsmaj.*2.*pi.*dogsmin);
                    
                    % Assert that the area of the 2D gaussian is the same
                    % select 2D gaussian crossing from center of RF
                    % assert(trapz(rf1(find(Y(:,1)==y0),:))==trapz(rf1(find(X(1,:)==x0),:)));
                    % assert(trapz(rf2(find(Y(:,1)==y0),:))==trapz(rf2(find(X(1,:)==x0),:)));
                    assert(isclose(trapz(rf1(find(Y(:,1)==y0),:)),trapz(rf2(find(Y(:,1)==y0),:)),'tolerance',0.0001),...
                                               'Area of gaussians not the same');
                    % Calculate the dog rf values
                    rf.values = rf1 - scale * rf2;
                    
                otherwise
                    error('%s not implemented yet', type)
            end
        end
        
        % Plot it
        function plot(rf,varargin)
            set(0, 'DefaultFigureRenderer', 'opengl');
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('rf'  ,  @(x)(isa(x,'pmRF')));
            p.addParameter('window',true, @islogical);
            p.parse(rf,varargin{:});
            w = p.Results.window;
            % Compute before plotting
            rf.compute
            % Plot it
            if w; mrvNewGraphWin('Receptive Field');end
            mesh(rf.values);
            grid on; 
            if w; xlabel('x'); ylabel('y'); end
            set(0, 'DefaultFigureRenderer', 'painters');
        end
    end
    
end



