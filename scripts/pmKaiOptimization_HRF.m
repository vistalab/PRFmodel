% Initial HRF constraints
% After this analysis see if we can contraint it more, try to find patterns on
%   the params

A1 = [5:1:10];
A2 = [5:1:10];
B1 = [.2:0.1:1];
B2 = [.2:0.1:1];
C  = [.1:0.1:.5];
minTotalWidth = 12;  % Total width from beginning to end of undershoot


% Create a grid of parameters and generate all possible HRFs
dontWork = [];
doWork   = [];
for a1=A1
    for a2=A2
        for b1=B1
            for b2=B2
                for c=C
                    % Generate time, TR=1sec
                    t   = 1:20;
                    % Generate HRF values (Friston two gammas)
                    hrf = (t/(a1*b1)).^a1   .* exp(-(t - (a1*b1))/b1) ...
                           - c*(t/(a2*b2)).^a2 .* exp(-(t-(a2*b2))/b2);
                                        
                    
                    % Add conditions:
                    
                    %  -1- Start with zero
                    startWithZero = isclose(hrf(1),0,'tolerance',0.01);
                    
                    %  -2- End with zero
                    endWithZero   = isclose(hrf(end),0,'tolerance',0.01);
                    
                    %  -3- There needs to be a undershoot
                    goesNegativeAt  = find(hrf < 0,1);
                    itHasUndershoot = ~isempty(goesNegativeAt);
                    
                    %  -4- The negative peak (undershoot) comes after the
                    %      positive peak
                    [~,maxPeakPos] = max(hrf);
                    [~,minPeakPos] = min(hrf);
                    peakOrder      = maxPeakPos < minPeakPos;
                    
                    %  -5- Overshoot needs to be wider than undershoot
                    %  -6- Total duration needs to be at least 12, after the
                    %      undershoot
                    if itHasUndershoot
                        % Calculate overshoot width
                        overshootWidth  = goesNegativeAt;
                        % Calcualte undershoot width
                        negativeHRF     = hrf(minPeakPos:end);
                        goesBackToZero  = find(negativeHRF > -0.01,1);
                        undershootWidth = (minPeakPos - goesNegativeAt) + goesBackToZero;
                        % Compare both
                        widerOvershoot  = overshootWidth > undershootWidth;
                        
                        % Calculate total length
                        hrfWidth        = minPeakPos + goesBackToZero;
                        itIsWideEnough  = hrfWidth >= minTotalWidth;
                    end
                    
                    % Modification to check a1b1<a2b2
                    justCheck = a1*b1 > a2*b2;
                    
                    
                    % - If it is a valid HRF, plot it
                    % - Store all valid combinations to extract patterns that can
                    %    be used as constraints
                    if startWithZero && endWithZero && itHasUndershoot && ...
                       peakOrder && widerOvershoot && itIsWideEnough && justCheck
                        plot(t,hrf,'color',[.6,.6,.6],'LineWidth',.5);hold on
                        doWork   = [doWork;[a1,a2,b1,b2,c]];
                    else
                        dontWork   = [dontWork;[a1,a2,b1,b2,c]];                        
                    end

                end
            end
        end
    end
end




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

