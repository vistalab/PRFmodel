
function st_makejson(prefix,savedir,temporalType,noiselevel,loadData,sigma,type,shuf)

%% load_default
defultJson = fullfile(stRootPath,'jsonfiles/default.json');
jsonText = fileread(defultJson);


%% set

if notDefined('prefix')
    prefix = [];
else
    prefix = strcat(prefix,'_');
end

if notDefined('savedir')
    savedir= './jsonfiles/';
else
    savedir = strcat(savedir,'/');
end


% loadData (load stimulus)
if notDefined('loadData')
    loadData = 'bar_9steps';
end

if notDefined('temporalType')
    temporalType = '2ch';
end

switch temporalType
    case 'div'
        temporalType = ['1ch-dcts'];
    case 'glm'
        temporalType = ['1ch-glm'];
    case '2ch'
        temporalType = ['2ch-exp-sig'];
end

if notDefined('noiselevel')
    noiselevel = 'low';
end

% which model the analysis is performed on
if notDefined('type')
    type ='st';
end

if notDefined('shuf')
    shuf=[0 1];
end

if notDefined('sigma')
    sigma=[1 6];
end



shuf = num2cell(shuf);
sigma = num2cell(sigma);

%%

if isfile(fullfile(pmRootPath,'data','stimulus',[loadData '.mat']))
    load(fullfile(pmRootPath,'data','stimulus',[loadData '.mat']))
    TR = num2str(size(stim,3));
else
   error("no stimfile found in the path") 
end

% TR = '120'; dur = '2';
% TR = '240'; dur = '5';

stimseq = {'a','b','c'}; % a b c
% stimseq = {'abc'}; % a b c

[stimseq, shuf,sigma]=BalanceFactors(1,0,stimseq,shuf,sigma);

%%

for ee =1 :length(stimseq)
%     jsonData = jsondecode(jsonText);
    J = loadjson(jsonText);
    if iscell(J)
        J=J{:};
    end
    jsonData = J;
    
    

    % stimulus
    if shuf{ee} == 0
        savename = [prefix type '_' TR stimseq{ee} num2str(sigma{ee}) '_' temporalType];
%         loadData = ['cst' '_' dur '_' TR];
        jsonData.Stimulus{1}.Shuffle = 0;
        jsonData.Stimulus{1}.myload = convertCharsToStrings(loadData);
    else
        savename = [prefix type '_' TR stimseq{ee} num2str(sigma{ee}) 's_' temporalType];
        stimName = [loadData '_shuffle'];
        jsonData.Stimulus{1}.Shuffle = 1;
        jsonData.Stimulus{1}.myload = convertCharsToStrings(stimName);

    end
    jsonData.Stimulus{1}.expName = convertCharsToStrings(savename);

    % overall
    jsonData.subjectName =  convertCharsToStrings(type);
    newStr = erase( convertCharsToStrings(savename) , convertCharsToStrings([type '_']));
    jsonData.sessionName = newStr;
    jsonData.Type = convertCharsToStrings(type);

    % sigma
    jsonData.RF{1}.sigmaMajor =sigma{ee};
    jsonData.RF{1}.sigmaMinor =sigma{ee};
        
    % noise
    if strcmp(noiselevel,'low')
        jsonData.repeats = 100;
        jsonData.Noise{1}.voxel = "low";
        jsonData.Noise{1}.seed = "random";
    elseif strcmp(noiselevel,'mid')
        jsonData.repeats = 100;
        jsonData.Noise{1}.voxel = "mid";
        jsonData.Noise{1}.seed = "random";
    elseif strcmp(noiselevel,'high')
        jsonData.repeats = 100;
        jsonData.Noise{1}.voxel = "high";
        jsonData.Noise{1}.seed = "random";
    elseif strcmp(noiselevel,'no')
        jsonData.repeats = 1;
        jsonData.Noise{1}.seed = "none";
    end
    
    % temporal
    jsonData.Temporal{1}.stimseq = convertCharsToStrings(stimseq{ee});
    jsonData.Temporal{1}.temporalModel = convertCharsToStrings(temporalType);
    
    % clearn up and change
    jsonString = jsonencode(jsonData);
    jsonString = strrep(jsonString, ',', sprintf(',\n'));
    jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
    jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));
    jsonString = [jsonString];
    fid = fopen([savedir savename '.json'], 'w');
    fwrite(fid, jsonString,'char');fclose(fid);
% 
    
end

%%
%     {'subjectName'      }
%     {'sessionName'      }
%     {'repeats'          }
%     {'useparallel'      }
%     {'TR'               }
%     {'Type'             }
%     {'signalPercentage' }
%     {'BOLDcontrast'     }
%     {'BOLDmeanValue'    }
%     {'computeSubclasses'}
%     {'HRF'              }
%     {'RF'               }
%     {'Stimulus'         }
%     {'Noise'            }



% jsonData.(fields{1}) = 'cst'