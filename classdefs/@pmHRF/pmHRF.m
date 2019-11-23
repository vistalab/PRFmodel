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
    % Plot several HRFs for comparison
    %{
      HRFs = {'vista_twogammas','popeye_twogammas','canonical','afni_spm'};
      pm = prfModel;
      pm.TR = 1.5;
      for HRF = HRFs
          pm.HRF.Type = HRF{:};
            
      end
    %}
    
    properties
        Type;          % Friston, Boyton, canonical...
        PM;            % prfModel that has some of the variables we need, such as TR
        Duration;      % Seconds
        params;        % Different values depending on the type of HRF
        values;
    end
    properties (Dependent = true, Access = public)
        tSteps;
    end
    properties(Dependent = true, SetAccess = private, GetAccess = public)
        % List of all possible options
        Types;         
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    %%
    methods (Static)
        function d = defaultsGet
            % Defaults for all different models
            d.Type                = 'friston';  % 'canonical'
            d.Duration            = 20;
            % Default params for friston
            % [peak1 width1 peak2 width2 -weight2]
            % [5.4000 5.2000 10.8000 7.3500 0.3500]
            % d.params.a            = [5.4  ,  5.2];
            % d.params.b            = [0.9, 0.9];
            % d.params.c            = 0.35;
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
        function v = get.Types(hrf)
            v = {'friston','boynton','canonical', 'vista_twogammas', ...
                         'popeye_twogammas', 'afni_gam','afni_spm'};
        end
        function compute(hrf)
            switch hrf.Type
                case {'random','rand','rnd'}
                    hrf.Type = hrf.Types{randi(length(hrf.Types))};
            end
            % Now with the new type, do the calculations. 
            switch hrf.Type
                case {'friston'}
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
                case {'boynton'}
                    % TODO: use a function under @HRF again, remove code from here.
                    params = hrf.params;
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
                    
                    t = hrf.tSteps;
                    
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
                case {'canonical'}
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
                    
                    % I am still not sure how it does the conversion to seconds.
                    % Inside analyzePRF it knows the TR, but for mrVista we need
                    % to provide the HRF in seconds.
                    % I had the code copied here, now use their function.
                    
                    
                    % Be careful here. getcanonicalhrf has the following params:
                    % thrf = getcanonicalhrf(hrf.params.stimDuration, hrf.TR);
                    % But when solving the 
                    thrf = getcanonicalhrf(hrf.TR, hrf.TR);
                    
                    
                    
                    % To be decided:
                    %   - do we want to be a fixed number of points?
                    %   - I can plot it correctly knowing the TR, but will it be
                    %     interpreted correctly by the different programs?
                    
                    % Add it to the output
                    hrf.values = thrf;
                case {'vista_twogammas'}
                    % This is the default inside the function
                    vistaParams     = [5.4 5.2 10.8 7.35 0.35];
                    hrf.values      = rmHrfTwogammas(hrf.tSteps, vistaParams);
                    pm = prfModel;
% plot examples
                    %                     pm.HRF.Type = 'vista_twogammas';
%                     figure(1);
%                     pm.TR = 1  ; subplot(3,1,1);pm.HRF.plot('window',false);
%                     pm.TR = 1.5; subplot(3,1,2);pm.HRF.plot('window',false);
%                     pm.TR = 2  ; subplot(3,1,3);pm.HRF.plot('window',false);
                case {'popeye_twogammas'}
                    % We obtain the values from python
                    
                    % I cannot make this work in GCP, so going back to hardocded
                    % Values copies from tests below
%                     
                    HRF1  = [0,    0.0153,    0.1804,    0.5041,    0.7814,    0.8771,    0.8018,    0.6336,    0.4445,    0.2745,    0.1371,    0.0320,   -0.0449,   -0.0977,   -0.1298,   -0.1440,   -0.1439,   -0.1335,   -0.1167,   -0.0970,   -0.0772,   -0.0591,   -0.0437,   -0.0314,   -0.0218,   -0.0148,   -0.0098,   -0.0064,   -0.0040,   -0.0025,   -0.0015,   -0.0009];
                    HRF15 = [0,    0.0706,    0.5041,    0.8541,    0.8018,    0.5384,    0.2745,    0.0808,   -0.0449,   -0.1162,   -0.1440,   -0.1397,   -0.1167,   -0.0870,   -0.0591,   -0.0372,   -0.0218,   -0.0121,   -0.0064,   -0.0032,   -0.0015,   -0.0007];
                    HRF2  = [0,    0.1804,    0.7814,    0.8018,    0.4445,    0.1371,   -0.0449,   -0.1298,   -0.1439,   -0.1167,   -0.0772,   -0.0437,   -0.0218,   -0.0098,   -0.0040,   -0.0015];
                    switch hrf.TR
                        case 1
                            hrf.values = HRF1;
                        case 1.5
                            hrf.values = HRF15;
                        case 2
                            hrf.values = HRF2;
                        otherwise
                            error('the HRF for this TR has not been generated')
                    end
                    
                    % hrf.values = double(py.popeye.utilities.double_gamma_hrf(0,hrf.TR));
                    % Check what is this function returning
                    %{
                    % I don't think this is right, I need to ask the developer
                    % if this is the expected behavior. 
                    % The shape of the hrf is very different just changing the
                    % TR
                    
                    HRF1  = double(py.popeye.utilities.double_gamma_hrf(0,1,1.0,''));
                    HRF15 = double(py.popeye.utilities.double_gamma_hrf(0,1.5,1.0,''));
                    HRF2  = double(py.popeye.utilities.double_gamma_hrf(0,2,1.0,''));
                    timeS1  = 0: 1  : 1  *(length(HRF1 )-1);
                    timeS15 = 0: 1.5: 1.5*(length(HRF15)-1);
                    timeS2  = 0: 2  : 2  *(length(HRF2 )-1);
                    figure(99); plot(timeS1,HRF1); hold on; plot(timeS15,HRF15);plot(timeS2,HRF2);
                    %}
                case {'afni_gam'}
                    hrfFileName = fullfile(pmRootPath,...
                        'data',['TR' num2str(hrf.TR) '_conv.ref.GAM.1D']);
                    if ~exist(hrfFileName,'file')
                        % Default GAM normalized to 1
                        system(['3dDeconvolve ' ...
                            '-nodata 50 ' num2str(hrf.TR) ' ' ...
                            '-polort -1 ' ... % Do not calculate detrending polinomials
                            '-num_stimts 1 ' ...
                            '-stim_times 1 "1D:0" GAM ' ... % k, tname, Rmodel
                            '-x1D ' hrfFileName]);
                    end
                    [~, hrf.values, ~, ~] = Read_1D (hrfFileName);
                case {'afni_spm'}
                    hrfFileName = fullfile(pmRootPath,...
                              'data',['TR' num2str(hrf.TR) 'sisar.conv.ref.SPMG1.1D']);
                    if ~exist(hrfFileName,'file')
                        system(['3dDeconvolve ' ...
                            '-nodata 50 ' num2str(hrf.TR) ' ' ...
                            '-polort -1 ' ...
                            '-num_stimts 1 ' ...
                            '-stim_times 1 "1D:0" SPMG1\(0\) ' ...
                            '-x1D ' hrfFileName]);
                    end
                    [~, hrf.values, ~, ~] = Read_1D (hrfFileName);
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
            switch hrf.Type
                case {'random','rand','rnd'}
                    hrf.Type = hrf.Types{randi(length(hrf.Types))};
            end
            % I think the biggest difference here is the interaction between
            % compute and tSteps:
            % --- some require to know tSteps to calculate hrf.values in compute
            % --- some require to know the number of values in hrf.values to be
            %     able to provide the correct tSteps. 
            % NOTE: if the Type is rand, it will set up one of the options
            % either the function get.tSteps or hrf.compute, the one that is
            % called first. 
            switch hrf.Type
                case {'1','friston','2','boynton','4','vista_twogammas'}
                    tSteps  = 0: hrf.TR: hrf.Duration;
                case {'3','canonical','5','afni_gam','6','afni_spm','7','popeye_twogammas'}
                    % We need to compute it first in order to know the length
                    hrf.compute;
                    tSteps  = 0: hrf.TR: hrf.TR*(length(hrf.values)-1);
                otherwise
                    error('HRF method %s not implemented or valid.',hrf.Type);
            end
            
        end
        
        % Plot it
        function plot(hrf,varargin)
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('hrf'  ,  @(x)(isa(x,'pmHRF')));
            p.addParameter('newwin', true, @islogical);
            p.parse(hrf,varargin{:});
            w  = p.Results.newwin;
            % Calculate it and return every time we need it.
            % Compute it just in case, to have the latest version
            hrf.compute;
            % Plot it
            if w; mrvNewGraphWin([hrf.Type ' HRF']);end;
            plot(hrf.tSteps, hrf.values,'-o');
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        end
        
    end
    
    
end



