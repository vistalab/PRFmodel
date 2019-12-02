function DTDT = pmForwardModelCalculate(DTDT)
% Creates a table with all parameters (defaults) required to perform a forward
% calculation
% 
%  Inputs: table with one or several rows of parameters to calculate bold series
%          in fwd model
% 
%  Outputs: returns the same table with the last column of pm-s calculated
% 
%  See also: forwardModelTableCreate
% 
%  GLU Vistalab 2019.05

% allpms = cell(height(DT),1);
% TODO: calculate chunksize dependnig on amount of cores parfor can take,
% with a maximum of 20000
chksize = 5000;

if height(DTDT) < chksize
    DTcc{1} = DTDT;  
    nchcks  = 1;
else
    nchcks = ceil(height(DTDT) / chksize);
    for nn=1:nchcks
        startindex = (nn*chksize) + 1 - chksize;
        if nn == nchcks
            endindex   = height(DTDT);
        else
            endindex   = nn*chksize;
        end
        DTcc{nn}   = DTDT(startindex:endindex,:);
    end
end
tic
parfor nn=1:nchcks
    DT = DTcc{nn};
    for ii=1:height(DT)
        % Do it row to row and parameter to parameter first, for debugging
        if mod(ii,100)==0
            disp([num2str(nn) ' -- ' num2str(ii) ' -- ' num2str(height(DT))])
        end
        %% Initialize the basic model with defaults
        dt     = DT(ii,:);
        % we need a fresh copy of the prfModel class here, otherwise it references
        % the same one and changes are not persistent
        pm     = prfModel;
        
        
        %% High Level Variables
        isprfmodel = @(x)(isa(x,'prfModel'));
        for vn = dt.Properties.VariableNames
            if ~istable(dt.(vn{:})) && ~isprfmodel(dt.(vn{:}))
                if iscell(dt.(vn{:}))
                    pm.(vn{:}) = dt.(vn{:}){:};
                else
                    pm.(vn{:}) = dt.(vn{:});
                end
            end
        end
        
        
        %% Stimulus
        for jj=1:width(dt.Stimulus)
            paramName               = dt.Stimulus.Properties.VariableNames{jj};
            pm.Stimulus.(paramName) = dt.Stimulus.(paramName);
        end
        % pm.Stimulus.compute;
        
        %% RF
        for jj=1:width(dt.RF)
            paramName         = dt.RF.Properties.VariableNames{jj};
            pm.RF.(paramName) = dt.RF.(paramName);
        end
        % pm.RF.compute;
        
        %% HRF
        for jj=1:width(dt.HRF)
            paramName          = dt.HRF.Properties.VariableNames{jj};
            val                = dt.HRF.(paramName);
            if iscell(val)
                pm.HRF.(paramName) = val{:};
            else
                pm.HRF.(paramName) = val;
            end
        end
        % pm.HRF.compute;
        
        
        %% Noise
        for jj=1:width(dt.Noise)
            paramName            = dt.Noise.Properties.VariableNames{jj};
            % voxel: we did not create it as a variable, as it is just to specify a
            % series of noise defaults. I don't like to have specific variable
            % names here, think about making it a noise param
            if strcmp(paramName,'seed')
                [val, status] = str2num(dt.Noise.(paramName){:});
                if status
                    pm.Noise.(paramName) = val;
                else
                    pm.Noise.(paramName) = dt.Noise.(paramName);
                end
            elseif ~strcmp(paramName,'voxel')
                pm.Noise.(paramName) = dt.Noise.(paramName);
            end
            
        end
        % pm.Noise.compute;
        
        %% Compute the synthetic signal
        % The compute at the top level computes all the lovel level ones.
        % Just do it once here.
        pm.compute;
        
        %% Assign it to the cell array (or Write back the updated pm model)
        % Parfor doesn't like this table
        % allpms{ii} = pm;
        DT.pm(ii) = pm;
        
    end
    DTcc{nn} = DT;
end
toc
%% Assign it back to the table before returning it.
% for ii=1:height(DT)
%     DT.pm(ii) = allpms{ii};
% end

%% Concatenate back
DTDT = DTcc{1};
if nchcks > 1
    for nn=2:nchcks
        DTDT = [DTDT; DTcc{1}];
    end
end



end
