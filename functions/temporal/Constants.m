classdef Constants
    % 	properties( Constant = true )
    %
    %     end
    %
    
    methods (Static)
        function d = getDir

            if isfile('/scratch/users/insubkim/yesSherlock.txt') 
                output_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/experiments/stRet/results/';
                simulator_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/experiments/simulator/';
            elseif isfile('/share/kalanit/users/insubkim/labComputer.txt')
                output_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/experiments/stRet/results/';
                simulator_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/experiments/simulator/';
            elseif isfile('/Volumes/Samsung_T5/yesSSD.txt')
                output_dir = '/Volumes/Samsung_T5/stRet/results';
                simulator_dir = '/Users/insubkim/oak_home/spatiotemporal/experiments/simulator/';
            else  
                output_dir = '/Users/insubkim/oak_home/spatiotemporal/experiments/stRet/results/';
                simulator_dir = '/Users/insubkim/oak_home/spatiotemporal/experiments/simulator/';
            end
            
 %%
              
 % test run --- test if the code runs or not
%  d.analysisoption = 99;
% %  d.analysisoption = 100;
%  sessionDate = 'test2';
 
 
% You can just disregard this one
            % sessionDate= '121420';

% first pass result, 2ch outperformed all while glm, DN showed similar results
% did it for subj01, and subj02 
%             sessionDate= '122320';    
%             d.analysisoption = 1;

% In a attempt to find out what is going on with the DN model
% I used different default params for the DN model
% from [0.1867 0 0.0895 3.5075 0.0647 0 1.0000] to  [0.05 0 0.1 2 0.1 0 1]
%             sessionDate = '010221';

% spatial smoothing
% analysisoption == 97
%             sessionDate = '011321'; % made a mistake? did not go through
%             d.analysisoption = 97;


% spatial smoothing 2
% Doing spatial smoothing again.....
% because I forgot to put smoothed estimated parameters back to the individual vertex
% now I need to plot it... but I've decided to work on the preds first            
%             d.analysisoption = 97;
%             sessionDate = '012021'; 
            
            
% I forgot to add "detrend" in the grid-fit so I fixed it
%             d.analysisoption = 97;
%             sessionDate = '020321'; 
% %             

            %  batch job 21005669
            % After bug fix re-running entirecode
%             d.analysisoption = 97;
%             sessionDate = '032421';
% 
% check after code cleaning if it works - > 
%             d.analysisoption = 99;
%             sessionDate = '032721';

%             
% % bug was not completly fixed... now lets check again
%             d.analysisoption = 97;
%             sessionDate = '032821';
%             
% now doing it with while saving each run preds
%             d.analysisoption = 97;
%             sessionDate = '040121';

            
%             bug was not completly fixed... now lets check again
%             d.analysisoption = 99;
%             sessionDate = 'test3';
% %             
%             d.analysisoption = 97;
%             sessionDate = 'simTest';
            
% run for all 3 subjects (subj01, 02, 03)
%             d.analysisoption = 97;
%             sessionDate = '032821';
%             
%             d.analysisoption = 99;
%             sessionDate = 'optim';

% 
%             % no longer doing coarse
%             d.analysisoption = 1;
%             sessionDate = '060921';

%             most updated Version
%             Using new ways of normalization

%             d.analysisoption = 1;
%             sessionDate = '070421';
% % 
%             d.analysisoption = 1;
%             sessionDate = '070721';
%             
            % version without crossvalidation
%             d.analysisoption = 3;
%             sessionDate = '081121';

%             d.analysisoption = 1;
%             sessionDate = '081221';

% %           this is with a normalize RF & linear based top 100
            d.analysisoption = 4;
            sessionDate = '081921';
            
            
%             d.analysisoption = 4;
%             sessionDate = '110821';

           

% 
%             d.analysisoption = 99;
%             sessionDate = 'test_070721';


            
            d.sessionDate = sessionDate;
            d.output_dir = output_dir;
            d.simulator_dir = simulator_dir;
            d.IRF_dir = fullfile(output_dir,'IRF',sessionDate,'/');
            d.grid_dir = fullfile(output_dir,'grid',sessionDate,'/');
            d.plot_dir = fullfile(output_dir,'plot',sessionDate,'/');
            
            
            % simulaiton related
            d.sim_json_dir = fullfile(simulator_dir,'jsonfiles',sessionDate,'/');
            d.sim_grid_dir = fullfile(simulator_dir,'grid',sessionDate,'/');
            d.sim_IRF_dir = fullfile(simulator_dir,'IRF',sessionDate,'/');
            d.sim_BOLD_dir = fullfile(simulator_dir,'synBOLD',sessionDate,'/');

            
        end

        function d = getTemporalParams
            
            d.tr = 1;     % 1 sec
            d.tnorm = 20; % transient channel normalization scalar
                       
            % zhou param, derived from the 2019 paper 
            % 'prm',   [0.1867 0 0.0895 3.5075 0.0647 0 1.0000] 
            % stig's zhou DN model defualt [0.05 0 0.1 2 0.1 0 1]
            
            d.temporalParams = {
            % 1ch-dcts
            % {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
            struct(...
            'type', '1ch-dcts', ...
            'fs', 1000, ... % sampling rate
            'fields',["tau1", "weight", "tau2", "nn", "delay", "shift", "scale"],...
            'num_channels', 1, ...
            'prm', [0.05 0 0.1 2 0.1 0 1], ...
            'searchRange',   [0.07  0  0.07 1  0.01 0 0; 1  0  1    6  0.5 0 0]...
            );  %   'prm', [0.05 0 0.1 2 0.1 0 1] ...

            
            %  2ch-exp-sig    
            %  {'tau_s', 'tau_ae', 'lambda_p', 'kappa_p', 'kappa_n', 'weight','shift'};
            %  searchRange: row1 = lowerRange, row2 = upperRange, 
            %               col  = params
            struct(...
            'type', '2ch-exp-sig', ...
            'fs', 1000, ...                         
            'num_channels', 2, ...,
            'fields', ["tau_s", "tau_ae", "Lp", "Kp", "Kn", "weight","shift"], ...
            'prm', [4.93 10000 0.1 3 3 0.5 0], ...
            'searchRange', [4  1 .01 .1 .1 0 0; 20 8 .5  7 7 0 0 ] ...    
            ); 
                    

            %  1ch-glm (linear)
            %  {'shift','scale'};
            struct(...
            'type', '1ch-glm', ...
            'fs', 1000, ...
            'fields',["shift", "scale"],...
            'num_channels', 1, ...
            'prm', [0 1], ...
            'searchRange',  [0 1; 0 1]...
            );
            };
        
        
        end
    end
    
    
    
end