classdef pmTemporal <   matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for Temporal attributes
    % Syntax:
    %
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
    
    %
    
    
    properties
        PM;         % prfModel that has some of the variables we need, such as TR
        temporalModel    ;  % [st] temporal modeling
        chan_preds       ;  % [st] saves channel-wise predictions
        run_preds        ;  % [st] saves channel-wise X RF predictions
        synBOLD          ;  % [st] saves final BOLD output with Noise added
        tParams          ;  % [st] temporal parameters
        fs               ;  % [st] temporal sample rate
    end
    
    properties (SetAccess = private, GetAccess = public)
        values;    % Result. Only can be changes from within this func.
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
        Name;
        nChan;
    end
    
    
    
    %%
    methods (Static)
        function d = defaultsGet
            d.temporalModel  = "None"; %stim temporal model
            d.fs = 1000;
            d.tParams.tParams = [];
            % Convert to table and return
            d = struct2table(d,'AsArray',true);
        end
    end
    methods
        % Constructor
        function temporal = pmTemporal(pm, varargin)
            % Obtain defaults table. If a parameters is not passed
            % it will use the default one defined in the static function
            d = temporal.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('pm',     @(x)(isa(x,'prfModel')));
            p.addParameter('temporalModel',d.temporalModel   , @ischar);
            p.addParameter('tParams',d.tParams   , @isstruct);
            p.addParameter('fs',d.tParams   , @isnumeric);
            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            temporal.PM             = pm;
            temporal.temporalModel    = p.Results.temporalModel;
            temporal.tParams          = p.Results.tParams;
        end
        
        function v = get.TR(temporal)
            v      = temporal.PM.TR;
        end
        
        function Name = get.Name(temporal)
            Name = [...
                'Exp-'          char(temporal.PM.Stimulus.expName), ...
                '_model-'       char(temporal.PM.Type),  ...
                '_tmodel-'      char(temporal.temporalModel) ...
                ];
        end
        
        function nChan = get.nChan(temporal)
            if contains(temporal.temporalModel,{'2ch','3ch','CST'})
                nChan = 2;
            else
                nChan = 1;
            end
        end        
        
        %% Methods available to this class and childrens, if any
        
        function compute(temporal)
            if strcmp(temporal.temporalModel,"None")
                % do nothing
            else
                temporal.PM.Stimulus.compute;
                temporal.PM.RF.compute;
                temporal.PM.HRF.compute;
                
                % get ms HRF
                hrf = resample(temporal.PM.HRF.values,temporal.fs,1)';
                hrf = hrf ./ sum(hrf(:));  % define TR
                tr = temporal.TR; % seconds
                
                % assign Temporal Model
                params.analysis.temporalModel = temporal.temporalModel;
                
                % getDefault temporal Param
                params = getTemporalParams(params);
                tParam = params.analysis.temporal.param;
                
                % over-ride temporal Param if userinput
                input_tParams= [];
                input_tParams = toSetField(input_tParams, cellstr(temporal.tParams.fields), temporal.tParams.values);
                if ~isempty(temporal.tParams)
                    MyFieldNames = fieldnames(tParam);
                    for fn = 1: length(MyFieldNames)
                        overlap = strcmp(fieldnames(input_tParams),MyFieldNames{fn});
                        if sum(overlap,'all')
                            tParam.(MyFieldNames{fn}) = input_tParams.(MyFieldNames{fn});
                        end
                    end
                    params.analysis.temporal.param = tParam;
                end
                
                params.analysis.hrf.func = hrf;
                params.analysis.temporal.param.fs = temporal.fs;
                params.analysis.temporal.fs = temporal.fs;
                params.analysis.spatial.values = temporal.PM.RF.values(:);
                params.verbose = 0;
                params.analysis.reluFlag = true;
                params.analysis.zeroPadPredNeuralFlag = true;
                params.analysis.normNeuralChan = false;
                params.analysis.normAcrossRuns = false;
                params.saveDataFlag = false;

                if strcmp(temporal.PM.Type,'st')
                    params.analysis.spatialModel = 'onegaussianFit';
                end
                
                [~,msStim] = st_keepWindowStim(temporal.PM,1);

                predictions = stPredictBOLDFromStim(params, msStim);
                
                temporal.chan_preds      = squeeze(predictions.predNeural);
                temporal.run_preds       = sum(predictions.predBOLD,3);

                
                
            end
        end
        
        
        
    end
end
