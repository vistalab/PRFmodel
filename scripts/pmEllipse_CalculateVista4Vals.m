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
tool        = 'vista4';
evc         = [bylabelsums.(tool).V1;bylabelsums.(tool).V2;bylabelsums.(tool).V3];

% CREATE ECCEN, ANGLE and FILTER 
[evc.angle,evc.eccen] = cart2pol(evc.x0, evc.y0);
evc.area              = pmEllipseArea(evc.sMaj,evc.sMaj);
% Look at the results and filter. 
%     SIZE : sMaj should be around 1<4
%     R2   : >25%
%     Eccen: 2<6
% Therefore, area will be: pi<50.27, because area=pi*sMin*sMaj
% But, looking at the distributions, it seems more reasonable to cut area at
% around 30
% Use this values when thresholding the mrVista HCP and synth

evc = evc(evc.r2>=0.25,:);
sMaj25 = prctile(evc.sMaj,25);
sMaj75 = prctile(evc.sMaj,75);
eccen25 = prctile(evc.eccen,25);
eccen75 = prctile(evc.eccen,75);

evc = evc(evc.sMaj>=sMaj25,:);
evc = evc(evc.sMaj<=sMaj75,:);
evc = evc(evc.eccen>=eccen25,:);
evc = evc(evc.eccen<=eccen75,:);

fprintf('\n Eccen at 25 = %.2g and 75 = %.2g percentiles.',eccen25,eccen75)
fprintf('\n Radius at 25 = %.2g and 75 = %.2g percentiles.',sMaj25,sMaj75)
fprintf('\n Area min = %.2g, max = %.2g\n',min(evc.area),max(evc.area))



kk = mrvNewGraphWin('vista4 values for synth select');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
nrow=2;ncol=2;
set(kk,'Position',[0.007 0.62  .9 .9]);
subplot(nrow,ncol,1)
plot(evc.x0,evc.y0,'b.')
axis equal
grid on
xlabel('X (deg)')
ylabel('Y (deg)')
title('pRF center positions (R2>25%)')

subplot(nrow,ncol,2)
histogram(evc.sMaj,25)
xlim([0,8])
xlabel('pRF Size (deg)')
grid on
title('pRF size histograms (R2>25%)')

subplot(nrow,ncol,3)
histogram(evc.eccen,25)
% xlim([-15,15])
grid on
xlabel('Eccentrcity (deg)')
title('pRF eccen histograms (R2>25%)')

subplot(nrow,ncol,4)
histogram(evc.area,25)
% xlim([-15,15])
xlabel('Area (deg^2)')
grid on
title('pRF area histograms (R2>25%)')

