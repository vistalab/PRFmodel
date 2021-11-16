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
    % DTcalc = synthBOLDgenerator(json, output_dir);
    
    % Examples
    %{
    %}
    
    %
    
    
    properties
        PM;         % prfModel that has some of the variables we need, such as TR
        %         stimseq          ;  % [st] exp A exp B  exp C
        temporalModel    ;  % [st] temporal modeling
        IRFpath          ;  % [st] path to save IRF predictions
        chan_preds       ;  % [st] saves channel-wise predictions
        run_preds        ;  % [st] saves channel-wise X RF predictions
        spc              ;  % [st] saves spc responses
        synBOLD          ;  % [st] saves final BOLD output with Noise added
        synBOLDpath      ;  % [st] path to save BOLD+noise predictions
        resolution       ;  % [st] Screen resolution (Hz)
        stim_on          ;  % [st] Number of on frames for each trial
        stim_off         ;  % [st] Number of off frames for each trial
        predur           ;
        tParams          ;  % [st] temporal parameters
        fs               ;
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
            d.resolution     = "60"; % Screen refresh Rate (Hz)
            d.predur = 0;
            d.fs = 1000;
            d.tParams.tParams = [];
%             d.tParams = st_temporal_default(d.temporalModel);
%             d.tParams = tparams.prm;

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
            p.addParameter('resolution',d.resolution   , @ischar);
            p.addParameter('predur',d.predur   , @double);
            p.addParameter('tParams',d.tParams   , @isstruct);

            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            temporal.PM             = pm;
            
            temporal.temporalModel    = p.Results.temporalModel;
            temporal.resolution       = p.Results.resolution;
            temporal.tParams          = p.Results.tParams;
            temporal.IRFpath          = Constants.getDir.sim_IRF_dir;
                        
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
           if contains(temporal.temporalModel,'2ch')
               nChan = 2;
           else
               nChan = 1;
           end
        end
        
%         function tparams = get.tParams(temporal)
%             tparams = st_temporal_default(temporal.temporalModel);
%             tparams = tparams.prm;
%         end

        
        %% Methods available to this class and childrens, if any
        
        function compute(temporal)
            if strcmp(temporal.temporalModel,"None")
                % do nothing
            else
                temporal.PM.Stimulus.compute;
                temporal.PM.RF.compute;
                temporal.PM.HRF.compute;
                hrf = resample(temporal.PM.HRF.values,1000,1)';

                
%                 c = Constants.getTemporalParams.temporalParams;
%                 
%                 for i = 1:length(c)
%                     if strcmp(c{i}.type, temp_type)
%                         idx = i;
%                     end
%                 end
                temp_type = char(temporal.temporalModel);
                temporal_param = temporal.tParams.tParams;
                if isempty(temporal_param)
                    tparams = st_temporal_default(temporal.temporalModel);
                    temporal_param = tparams.prm;
                end
                num_channels   = temporal.nChan;
                
                % define ms hrf
%                 hrf = canonical_hrf(1 / temporal.fs, [5 14 28]);

                % define TR
                tr = Constants.getTemporalParams.tr; % seconds
                
                % compute rfStim
                [keep,msStim] = st_keepWindowStim(temporal.PM,0);
                rf = temporal.PM.RF.values(keep);
                rfStim = full(msStim'*sparse(rf));
                
                % account and clip the preduration period
%                 rfStim = rfStim(temporal.predur*1000+1:end);

                
                t = 0.001 : 0.001 : size(rfStim,1)/temporal.fs;
                
                % compute rsp
                rsp =st_tModel(temp_type,temporal_param,rfStim', t); % rsp = time (ms) x space
                
                % account and clip the preduration period
                %                 rsp = cellfun(@(x) x(temporal.predur*1000+1:end,:),rsp,'UniformOutput',false);
                
                
                temporal.run_preds = cellfun(@(X, Y) convolve_vecs(X, Y, temporal.fs, 1 /tr), ...
                    rsp, repmat({hrf}, size(rsp)), 'uni', false);
                temporal.run_preds = cell2mat(temporal.run_preds);
                temporal.run_preds = temporal.run_preds(:,1:num_channels);
                
                %                 normMax = @(x) x./max(x);
                %                 run_preds = normMax(run_preds);
                %
                weight = 1; 
                if contains(temp_type,'2ch')
                    maxNorm = max(temporal.run_preds(:,1))/ max(temporal.run_preds(:,2));
                    temporal.run_preds(:,2) = temporal.run_preds(:,2)*maxNorm;
                    weight = 0.5;
                end
                
                temporal.chan_preds      = rsp;
                temporal.run_preds       = temporal.run_preds./max(temporal.run_preds)*weight;

%                 temporal.run_preds       = temporal.run_preds*weight;
            end
        end
        
        
        
    end
end



