function pmEllipse_CalculateVista4Vals
% Make Figures 6A-B-C, S7 and S8
% 
% TODO: separate it in sub-scripts that use the same dataset that we load here
% 
% See also
%  s00_MainFiguresScript
% 
%% Plotting parameters
clear all
ext  = 'png'; % Could be svg
saveTo = fullfile(pmRootPath,'local','figures');  % Folder path
if ~exist(saveTo,'dir'), mkdir(saveTo); end


%% READ: Real Data 7T

proj   = 'realdata';
tools  = {'vista4'};
subs   = {'115017','164131','536647'}; 
ses    = '01';
run    = '01';

[compTable,bylabelsums] = pmEllipse_loadExpData(proj,tools,subs,ses,run);
nonfilteredbuylabelsums = bylabelsums;

%% Calculate the values
bylabelsums = nonfilteredbuylabelsums;
tool = 'vista4';
evc  = [bylabelsums.(tool).V1;bylabelsums.(tool).V2;bylabelsums.(tool).V3];

kk = mrvNewGraphWin('vista4 values for synth select');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5 0.4]);
subplot(1,2,1)
plot(evc.x0,evc.y0,'b.')
axis equal
grid on
title('pRF center positions')

subplot(1,2,2)
histogram(evc.sMaj,100)
xlim([0,8])
grid on
title('pRF size histograms')

kk = mrvNewGraphWin('vista4 values for synth select');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5 0.4]);
subplot(1,2,1)
histogram(evc.x0,100)
xlim([-15,15])
grid on
title('pRF x0 histograms')

subplot(1,2,2)
histogram(evc.y0,100)
xlim([-15,15])
grid on
title('pRF y0 histograms')

