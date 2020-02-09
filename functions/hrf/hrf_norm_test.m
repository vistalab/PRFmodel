

doIt     = false;
window   = false;
plotbold = true;
plotstim = false;
plotrf   = false;
ylims    = [-1,2];
xlims    = [0,70];
viewn    = 2;
tr       = 1.5;
hrtype   = 'boynton';
hrfnorm  = 'norm';
spc      = 'frac';

hh = mrvNewGraphWin('HRF comparison');
set(hh,'Position',[0.007 0.62  0.4  0.4]);
nrows = 1; ncols = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the middle bar 20 frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(nrows,ncols,1)
pm = prfModel;
pm.TR=tr;
pm.RF.sigmaMajor = .25;
pm.RF.sigmaMinor = .25;
pm.HRF.Type      = hrtype;
pm.HRF.normalize = hrfnorm;
pm.signalPercentage=spc;
pm.Stimulus.durationSecs = 40;
pm.Stimulus.compute
if plotstim; pm.Stimulus.plot('window',window); end
if doIt
    % Edit it manually to make it a single on/off bar
    fName = fullfile(pmRootPath,'data','stimulus',['Exp-103_bin-true_size-20x20_resize-true_Horz-101x101_barW-2_dur-40_TR-' num2str(tr) '_framedur-4_Shuffle-false.mat']);
    kk = load(fName);
    % Repeat fourth bar, make rest zero, and save it with the same name, it will
    % not overwrite it.
    stim = kk.stim;
    stim(:,:,3:22) = repmat(kk.stim(:,:,5),[1,1,20]);
    stim(:,:,[1,2,23:end]) = zeros(size(kk.stim(:,:,[1,2,23:end])));
    % Check it
    montage(stim)
    % Save it
    save(fName,'stim');
    % Check again as part of class, it should just read it and plot what we want
    pm.Stimulus.plot
end
% Compute the time series
% pm.compute: we can't do this because the number of basis functions too low, for slow drift noise
pm.computeBOLD
if plotbold
    pm.plot('what','nonoisetimeseries','window',window)
    xlim(xlims)
    ylim(ylims)
    text(0,-0.5,sprintf(' No scaled BOLD, sum: %0.2f, mean: %0.2f, pow: %0.2f, contrast: %0.2f', ...
                         sum(pm.BOLDconv), ...
                         mean(pm.BOLDconv), ...
                         rms(pm.BOLDconv)^2, ...
                         (max(pm.BOLDconv)-min(pm.BOLDconv))/2), ...
         'FontSize',14)
     text(0,-0.7,sprintf(' Scaled BOLD, sum: %0.2f, mean: %0.2f, pow: %0.2f, contrast: %0.2f', ...
                         sum(pm.BOLD), ...
                         mean(pm.BOLD), ...
                         rms(pm.BOLD)^2, ...
                         (max(pm.BOLD)-min(pm.BOLD))/2), ...
         'FontSize',14)
end
if plotrf;pm.RF.plot('window',window);view(viewn);axis equal;
text(0,-2,sprintf('Sum of RF values: %.2f, vol:%.1f',sum(pm.RF.values(:)),...
    trapz(pm.Stimulus.XY{2}(:,1),trapz(pm.Stimulus.XY{1}(1,:),pm.RF.values,2),1)))
end
% title(['middle bar 20 frames, RF .25 deg center, TR=' num2str(tr)])
title(['HRF type: ' hrtype ', norm: ' hrfnorm])
