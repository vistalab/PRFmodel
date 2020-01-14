%{
clear all; close all; clc
p = '/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL';
f = 'sub-ellipse_ses-sess02-prf_acq-normal_run-01_bold.mat';
load(fullfile(p,f))
%}



kk = mrvNewGraphWin('NoiselessCloudPoints','wide');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.2  0.3]);
% subplot(1,4,1)
tools  = {'afni'};
useHRF = 'afni_spm';
nslvl  = 'mid';
pmCloudOfResults(compTable   , tools ,'onlyCenters', false ,'userfsize' , 2, 'centerdistr',false,...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'useellipse',true, 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',true, ...
                 'color', [0.5,0.5,0.5], 'xlims',[-3, 6],'ylims',[-3,6],...
                 'addcihist', false, 'xtick',[-2:6],'ytick',[-2:6],...
                 'location', [3,0], ... % 'all', ... % [3,3], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
             
             
             
%% NOTES ON THETA
% AFNI
% This is how we read it:
pmEstimates.Theta = results(:,6);
% This is the explanation of Theta by AFNI
%                   Given stimulus images over time s(x,y,t), find x0, y0, sigma, R and
%                   theta values that produce a best fit of the model to the
%                   data.  Here x0, y0 are taken to be the center of the
%                   population receptive field, sigma is the minor width of it
%                   (sigma_x, below), sigrat R is the ratio (sigma_y / sigma_x),
%                   and theta is the rotation from the y-direction major axis
%                   (so zero is in the positive y-direction).
% 
%                   We assume sigma_y >= sigma_x and refer to sigrat >= 1, since that
%                   sufficiently represents all possibilities.  The reciprocol would
%                   come from the negative complimentary angle, and would therefore be a
%                   redundant solution.
% 
%                    parameter domains:
%                      x,y        : [-1,1], scaled by the mask, itself
%                      sigma      : (0,1], where 1 means the mask radius
%                      R (sigrat) : [1,inf), since sigma defines the smaller size
%                      theta      : [-PI/2, PI/2), since rotation by PI has no effect





























