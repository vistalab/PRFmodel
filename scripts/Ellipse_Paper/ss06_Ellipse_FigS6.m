function ss06_Ellipse_FigS6(saveTo,ext)
if ~isfolder(saveFigTo); mkdir(saveFigTo); end


% NoiseCircleTest (Simulation:  Take a circle and make two noisy estimates of
% its diameter.  (D1 + noise)/(D2 + noise).  This ratio will always be centered
% on 1.  But now, use the same data to estimate
% (max(D1+noise,D2+noise))/min((D1+noise,D2+noise).  This will always be > 1.
% How much greater?  If your estimate is < (1.2?  1.3?) then the data are
% consistent with a circle.

% Noisy Circle Tests.m
% We have some results in mrVista with simulated data that are not perfect but
% that can be explained by the very nature of the calculation. 

noiseFactor  = .7; % Std. Dev.
radiuses     = [0.25:0.25:6]';  % degrees
% freeRatiosMn = zeros(size(radiuses)); 
% freeRatiosSt = zeros(size(radiuses)); 
limRatiosMed  = zeros(size(radiuses)); 
limRatiosMin  = zeros(size(radiuses)); 
limRatiosMax  = zeros(size(radiuses)); 
for nr=1:length(radiuses)
    radius = radiuses(nr);
    % Generate the random two radiuses
    rng(44444,'twister')
    n1    = noiseFactor * randn(1000,1);
    R1    = radius + n1;
    R1ind = R1 > 0.01;
    
    rng(54321,'twister')
    n2    = noiseFactor * randn(1000,1);
    R2    = radius + n2;
    R2ind = R2 > 0.01;
    
    % Combine all positive radius (AND)
    Rind = R1ind & R2ind;
    R1   = R1(Rind);
    R2   = R2(Rind);

    % Calculate the ratio that is always positive
    limRatio         = max([R1,R2],[],2) ./ min([R1,R2],[],2);
    
    % Calculate median, min and max
    limRatiosMed(nr)  = median(limRatio);
    limRatiosMin(nr)  = min(limRatio);
    limRatiosMax(nr)  = max(limRatio);
end
    
% ext       = 'svg';
fnameRoot = 'FigS6_RadiusSims';
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.4  0.4]);

ystart =     ones(size(radiuses));
ystop  = 6 * ones(size(radiuses));
plot([radiuses.';radiuses.'],[ystart.';ystop.'], ...
    'LineWidth',.7,'LineStyle','-.','Color','k')
hold on
% plot([0,6],[1,1],'LineWidth',1.5,'LineStyle','--','Color','k') % 0.75*[0 1 0])
% Cs              = 0.65*distinguishable_colors(1+length(sizes),'w');
% Apply percentiles and plot individually
for ne=1:length(radiuses)
    rad          = radiuses(ne);
    limRatioMed  = limRatiosMed(ne);
    limRatioMin  = limRatiosMin(ne);
    limRatioMax  = limRatiosMax(ne);
        
    % Plot it
    b   = scatter(rad,limRatioMed,80,'ko','filled');
    vax = plot(rad * [1,1],[limRatioMin, limRatioMax], ...
             'Color','k','LineStyle',':','LineWidth',3);    
end

% legend([a,b], {'Free Ratios','Limited Ratios'})
title(strrep(fnameRoot,'_','\_'))
xlabel('Radius')
ylabel('Aspect Ratio')
ylim([0,8]);
set(gca, 'FontSize', 16)
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);



end

