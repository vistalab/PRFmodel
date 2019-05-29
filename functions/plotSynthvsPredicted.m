function gcf = plotSynthvsPredicted(pm, result, varargin)
% Takes a pm and its corresponding prf model fit result and plots them
% 
%  Inputs: pm of the synthetic data and the result
%           TODO: right now only does analyzePRF
% 
%  Outputs: gcf
% 
% 

%% Read parameters
varargin = mrvParamFormat(varargin);
            
p = inputParser;
p.addRequired('pm');
p.addRequired('result');
p.addParameter('saveplot'     , true         , @islogical); % To save images in plotfilename
p.addParameter('plotfilename' , './plot.png' , @ischar);  % plotfilename

p.parse(pm,result,varargin{:});
            
saveplot     = p.Results.saveplot;
plotfilename = p.Results.plotfilename;

%% Create the plot

stimulus = {stimValuesRead(pm.stimulus.values)};
data     = {pm.BOLD.predictedWithNoise};

% Define some variables
res   = [pm.stimulus.fieldofviewVert, pm.stimulus.fieldofviewHorz]; % row x column resolution of the stimuli
resmx = max(res);                   % maximum resolution (along any dimension)
hrf   = result.options.hrf;         % HRF that was used in the model
degs  = result.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Prepare the stimuli for use in the model
stimulusPP = {};
for p=1:length(stimulus)
  stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

% Define the model function.  This function takes parameters and stimuli as input and
% returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
% of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
% Although it looks complex, what the function does is pretty straightforward: construct a
% 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
% Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
% taking care to not bleed over run boundaries.
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

% Construct projection matrices that fit and remove the polynomials.
% Note that a separate projection matrix is constructed for each run.
polymatrix = {};
for p=1:length(degs)
      polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
end


% Which voxel should we inspect?  Let's inspect the second voxel.
vx = 1;  % By definition, we will have only one (each row)

% For each run, collect the data and the model fit.  We project out polynomials
% from both the data and the model fit.  This deals with the problem of
% slow trends in the data.
datats = {};
modelts = {};
for p=1:length(data)
  datats{p} =  polymatrix{p} * data{p}(vx,:)';
  modelts{p} = polymatrix{p} * modelfun(result.params(1,:,vx),stimulusPP{p});
end

% Visualize the results
gcf = mrvNewGraphWin('Synthetic vs AnalyzePRF'); hold on;
set(gcf,'Units','points','Position',[100 100 1000 400]);
plot(cat(1,datats{:}),'r-'); hold on;
plot(cat(1,modelts{:}),'b-');
% straightline(300*(1:4)+.5,'v','g-');
xlabel('Time (s)');
ylabel('BOLD signal');
ax = axis;
% axis([.5 1200+.5 ax(3:4)]);
title('Time-series data');


%% Save the plot
    if saveplot
        [FILEPATH,NAME,EXT] = fileparts(plotfilename);  % png or svg
        saveas(gcf,plotfilename, EXT);
    end
end

