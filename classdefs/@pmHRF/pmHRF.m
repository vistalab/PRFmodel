classdef pmHRF <  matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for HRF-s. Every particular instance of this class
    % will have different parameters, so it will be a children class. For
    % example, "Friston" or "Boynton"
    %
    % Syntax:
    %      hrf = HRF;
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
        Type;          % Firston, Boyton, canonical...
        PM;            % prfModel that has some of the variables we need, such as TR
        Duration;      % Seconds
        params;        % Different values depending on the type of HRF
        values;
    end
    
    properties (Dependent = true, Access = public)
        tSteps;
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    %%
    methods (Static)
        function d = defaultsGet
            % Defaults for all different models
            d.Type                = 'friston';
            d.Duration            = 20;
            % Default params for friston
            d.params.a            = [6  ,  12];
            d.params.b            = [0.9, 0.9];
            d.params.c            = 0.35;
            % Default params for friston
            d.params.n            = 3;
            d.params.tau          = 1.08;
            d.params.delay        = 2.05;
            % Default params for canonical
            d.params.stimDuration = 1;
            
            % Convert to table and return
            d = struct2table(d,'AsArray',true);
        end
    end
    methods
        function hrf = pmHRF(pm,varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = hrf.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired('pm'       ,             @(x)(isa(x,'prfModel')));
            p.addParameter('type'    ,d.Type{:}  , @ischar);
            p.addParameter('params'  ,d.params   , @isstruct);
            p.addParameter('duration',d.Duration , @isnumeric);
            p.parse(pm,varargin{:});
            % Assign it
            params   = p.Results.params;
            Duration = p.Results.duration;
            
            % Initialize the pm model and hrf model parameters
            hrf.PM       = pm;
            hrf.Type     = p.Results.type;
            hrf.Duration = p.Results.duration;
            % Check if only some of the fields in params where set
            hrf.params   = pmParamsCompletenessCheck(p.Results.params,d.params);
        end
        
        function compute(hrf)
            switch hrf.Type
                case 'friston'
                    % TODO: use a function under @HRF again, remove code from here.
                    a = hrf.params.a;
                    b = hrf.params.b;
                    c = hrf.params.c;
                    
                    t = hrf.tSteps;
                    % Calculate d
                    for ii = 1:2
                        d(ii) = a(ii)*b(ii);
                    end
                    % Calculate actual values
                    hrf.values = (t/d(1)).^a(1)   .* exp(-(t - d(1))/b(1)) ...
                        - c*(t/d(2)).^a(2) .* exp(-(t-d(2))/b(2));
                case 'boynton'
                    % TODO: use a function under @HRF again, remove code from here.
                    
                    % Make sure that when boynton is selected, params comes with
                    % the correct parameters, otherwise use defaults.
                    if ~isfield(params,'n'),    params.n     = 3;   end
                    if ~isfield(params,'tau'),  params.tau   = 1.08;end
                    if ~isfield(params,'delay'),params.delay = 2.05;end
                    % Update the obj just in case
                    hrf.params = params;
                    % Assign them internally to be used in the function
                    n     = hrf.params.n;
                    tau   = hrf.params.tau;
                    delay = hrf.params.delay;
                    
                    %    hrf = boyntonHIRF(t, [n=3], [tau=1.08], [delay=2.05]])
                    %
                    %Purpose:
                    %   Compute the Boynton et al. HIRF function.  The return values include
                    % both the hrf and the values of time and the parameters used in the
                    % computation.
                    %
                    % Equation:  (Eq 3 from Boynton & Heeger, J Neurosci 1996)
                    %   h(t) = [(t/tau) ^ (n-1) * exp(-t/tau)] / [tau(n-1)!]
                    %
                    % Inputs:
                    %   t: t should be in units of SECONDS, and reflect the time window in
                    %   which to estimate the HRF.
                    %
                    %   n: exponent in eq. above. (corresponds to z in standard gamma
                    %   functions).
                    %
                    %   tau: time constant in eq.
                    %
                    %   delay: additional delay before onset of gamma function. This is an
                    %   added heuristic, which may vary from subject to subject.
                    %
                    %
                    % Outputs:
                    %   hrf: estimate of the HRF sampled at t seconds.
                    %   IMPORTANT: If you are going to convolve
                    %   this for fMRI analysis, and the TR of your data is not 1, you will
                    %   need to resample both the input t and hrf to match the MR frames.
                    %
                    % written 2005 by wandell.
                    % ras, 01/2007: heavily modified: no parms struct, each arg can be
                    % specified separately; clarifies difference between seconds and frames,
                    % doesn't modify t.
                    % disp('Boynton HIRF')
                    
                    
                    % initialize the HRF to be zeros, the same size as t:
                    thrf = zeros( size(t) );
                    
                    % The HRF is not specified for t < 0 secs. In addition, any
                    % values before the delay should also be zero.
                    % So, we only sample the HRF below, for time points after [0+delay].
                    % We call this sampling vector x to distinguish it from t.
                    x = t - delay;
                    x = x(x>0);  % not defined for x<0
                    
                    if round(n)~=n
                        warning('Boynton HIRF parameters stored as unusual values.  n is being rounded from %.2f to %i.',n, round(n))
                        disp('Check stored glmHRF_params against cannonical values.')
                        n = round(n);
                    end
                    
                    % main computation (per Boynton & Heeger, 1996)
                    tmpHrf = (x/tau).^(n-1) .* exp(-(x/tau)) / (tau*(factorial(n-1)));
                    
                    % paste in tmpHRF into appropriate indices in hrf, corresponding
                    % to the (shifted after delay) time points:
                    thrf(end-length(tmpHrf)+1:end) = tmpHrf;
                    % Assign it to the output of our object
                    hrf.values = thrf;
                case 'canonical'
                    % average empirical HRFs from various datasets (where the HRFs were
                    % in response to a 3-s stimulus) and then do the requisite
                    % deconvolution and convolution to obtain a predicted HRF to a
                    % stimulus of duration <duration>, with data sampled at a TR of <tr>.
                    %
                    % the resulting HRF is a row vector whose first point is
                    % coincident with stimulus onset.  the HRF is normalized such
                    % that the maximum value is one.  note that if <duration> is
                    % small, the resulting HRF will be quite noisy.
                    %
                    % example:
                    % hrf = getcanonicalhrf(4,1);
                    % figure; plot(0:length(hrf)-1,hrf,'ro-');
                    
                    % load HRFs from five datasets and then take the average.
                    % these were the empirical response to a 3-s stimulus, TR 1.323751 s
                    % hrf = mean(catcell(2,getsamplehrf([9 10 11 12 14],1)),2)';  % 1 x time
                    % store a hard copy for speed
                    thrf    = [0 0.0314738742235483 0.132892311247317 0.312329209862644 0.441154423620173 0.506326320948033 0.465005683404153 0.339291735120426 0.189653785392583 0.0887497190889423 0.0269546540274463 -0.00399259325523179 -0.024627314416849 -0.0476309054781231 -0.0550487226952204 -0.0533213710918957 -0.0543354934559645 -0.053251015547776 -0.0504861257190311 -0.0523878049128595 -0.0480250705100501 -0.0413692129609857 -0.0386230204112975 -0.0309582779400724 -0.0293100898508089 -0.0267610584328128 -0.0231531738458546 -0.0248940860170463 -0.0256090744971939 -0.0245258893783331 -0.0221593630969677 -0.0188920336851537 -0.0205456587473883 -0.0230804062250214 -0.0255724832493459 -0.0200646133809936 -0.0101145804661655 -0.014559191655812];
                    trorig = 1.323751;
                    
                    
                    % Obtain the parameters for this instance
                    duration = hrf.params.stimDuration;
                    tr       = hrf.TR;
                    
                    
                    % resample to 0.1-s resolution
                    % GLU: added the condition for 1.82. Generalize it for other cases.
                    trnew = 0.1;
                    if tr==1.82
                        trnew=0.02;
                    end
                    thrf = interp1((0:length(thrf)-1)*trorig,thrf,0:trnew:(length(thrf)-1)*trorig,'PCHIP');
                    
                    % deconvolve to get the predicted response to 0.1-s stimulus
                    thrf = deconvolvevectors(thrf,ones(1,3/trnew));
                    
                    % convolve to get the predicted response to the desired stimulus duration
                    thrf = conv(thrf,ones(1,duration/trnew));
                    
                    % resample to desired TR
                    thrf = interp1((0:length(thrf)-1)*trnew,thrf,0:tr:(length(thrf)-1)*trnew,'PCHIP');
                    
                    % make the peak equal to one
                    thrf = thrf / max(thrf);
                    plot(thrf)
                    % make the same time points as TR
                    thrf = thrf(1:length(hrf.tSteps));
                    
                    
                    % Add it to the output
                    hrf.values = thrf;
                otherwise
                    error('HRF method %s not implemented or valid.',hrf.Type);
            end
        end
        
        
        % Methods available to this class and his childrens (friston, boynton... classes)
        function v = get.TR(hrf)
            v = hrf.PM.TR;
        end
        function tSteps = get.tSteps(hrf)
            % Calculate it and return every time we need it.
            % TODO: Does it need to be in TR intervals?
            tSteps  = 0: hrf.TR: hrf.Duration;
        end
        
        
        
        
        % Plot it
        function plot(hrf)
            % Calculate it and return every time we need it.
            % Compute it just in case, to have the latest version
            hrf.compute;
            % Plot it
            mrvNewGraphWin([hrf.Type ' HRF']);
            plot(hrf.tSteps, hrf.values);
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        end
        
    end
    
    
end



