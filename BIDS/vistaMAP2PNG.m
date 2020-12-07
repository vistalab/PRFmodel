function vistaMAP2PNG(projectDir, subject,session,runnums)
% Visualize an MGZ ret data on a freesurfer surface in Matlab
%
% Maps2PNG(bidsfolder, subject, session, desc)
%
% INPUTS
%   bidsfolder     : path to BIDS project
%   subject        : BIDS subject name
%   session        : BIDS session name
%   desc           : name of subfolder reflecting model/design that was
%                    used for analysis;

% Example
% bidsfolder = '/Volumes/server/Projects/SampleData/BIDS/';
% subject    = 'wlsubj042';
% session    = '01';
% desc       = 'coarse';



% set up our paths

% latestDir = find_latest_dir(projectDir);
latestDir  = projectDir;
resultsdir = fullfile (projectDir,'derivatives','prfanalyze-vista06', ...
             sprintf('sub-%s',subject), sprintf('ses-%s',session));
hcpdir     = '/data/localhome/glerma/toolboxes/PRFmodel/local/realdata/hcp';

% fspth = fullfile(projectDir, 'derivatives', 'freesurfer',  ['sub-' subject], 'surf');
fspth      = fullfile(hcpdir, 'fsaverage', 'surf');
figureDir  = fullfile(projectDir, 'figures');


% if the fig dir does not exist, make one
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

import mlreportgen.report.*
import mlreportgen.dom.*

rpt = Report(sprintf('%s/pRFSummary_%s',figureDir,subject));
p1  = Paragraph(sprintf('prFSummary'))
p1.Bold=true;
p1.FontSize = '30';
add(rpt,p1);


p1 = Paragraph(sprintf('This is a report for %s with %.0f repetitions',subject,length(runnums)));
p1.FontSize = '20';
add(rpt,p1);

hemispheres = {'lh';'rh'};
path2roi = {'V1_exvivo';'V2_exvivo'};

% if the fig dir does not exist, make one
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end


% find the prf data, assuming we want images of everything
mapsList = {'angle', 'eccen', 'sigma', 'vexpl'};
map_file = [];

% load all the data

for hemi = 1 : length(hemispheres)
    
    for thisMap = 1:length(mapsList)
        
        map_file.(mapsList{thisMap}).(hemispheres{hemi}) = load_mgh(fullfile(resultsdir, sprintf('%s.%s.',hemispheres{hemi},mapsList{thisMap})));
        
    end
    
end

% loop through the maps and create png pics

chapter = Chapter();
chapter.Title = 'pRFMaps';
add(rpt,chapter);

for thisFig = 1:length(mapsList)
    
%     fig=figure(1);clf
    fig = figure(1);clf
    
    for hemi = 1 : length(hemispheres)
        
        roi = [];
        
        for r = 1 : length(path2roi)
            
            ind  = read_label(['sub-' subject],sprintf ('%s.%s%s',hemispheres{hemi},path2roi{r}));
            roi  = [roi; ind(:,1) + 1];
            
        end
        
        curv = read_curv(sprintf('%s/%s.curv',fspth,hemispheres{hemi}));
        surf_file = fullfile(fspth, sprintf('%s.inflated',hemispheres{hemi}));
        [vertices, faces] = read_surf(surf_file);
        
        myroi = zeros(size(map_file.vexpl.(hemispheres{hemi})));
        myroi(roi) = 1;
        
        
        thr  = double(map_file.vexpl.(hemispheres{hemi})>0) & double(map_file.eccen.(hemispheres{hemi})<Inf) & myroi ;
        
        subplot(2,2,hemi)
        plot_mesh(faces, vertices, curv, 'gray');
        freezeColors
        mymap = map_file.(mapsList{thisFig}).(hemispheres{hemi}).*thr;
        mymap(mymap==0) = NaN;
        
        variable = mapsList{thisFig};
        
        switch variable
            
            case 'angle'
                
                plot_mesh(faces, vertices, mymap,'hsv');
                caxis([-pi pi])
                
            case 'eccen'
                
                plot_mesh(faces, vertices, mymap, 'jet');
                caxis([0 12])
                
            case 'sigma'
                
                plot_mesh(faces, vertices, mymap, 'parula');
                caxis([0 3])
                
                
            case 'vexpl'
                
                plot_mesh(faces, vertices, mymap, 'hot');
                caxis([0 1])
                
        end
        
        if hemi == 1
            set_view(gcf)
            view(45,0)
            
        else
            set_view(gcf)
            view(-45,0)
            
        end
        
    end
    
    subplot(2,2,3:4)
    
    switch variable
        
        case 'angle'
            
            cbh=colorbar('North');
            colormap(hsv)
            caxis([-pi pi])
            title(variable)
            set(cbh,'XTick',[-pi 0 pi])
            set(cbh,'XTickLabel',{'-\pi';'0';'\pi'})
            
        case 'eccen'
            
            cbh=colorbar('North');
            colormap(jet)
            caxis([0 12])
            title(variable)
            set(cbh,'XTick',[0 4 8 12])
            set(cbh,'XTickLabel',{'0';'4';'8';'12'})
            
        case 'sigma'
            
            cbh=colorbar('North');
            colormap(parula)
            caxis([0 3])
            title(variable)
            set(cbh,'XTick',[0 1 2 3])
            set(cbh,'XTickLabel',{'0';'1';'2';'3'})
            
        case 'vexpl'
            
            cbh=colorbar('North');
            colormap(hot)
            caxis([0 1])
            title(variable)
            set(cbh,'XTick',[0 0.25 0.5 0.75 1])
            set(cbh,'XTickLabel',{'0';'0.25';'0.5';'0.75';'1'})
            
            
    end
    
    
  
%     fig.Snapshot.Caption = mapsList{thisFig};
%     fig.Snapshot.Height = '5in';

    axis off
    set(gca,'Fontsize',20)
    
    figReporter0 = Figure(fig);
%     add(rpt,figReporter0);
%     t_map{thisFig} = Image(getSnapshotImage(figReporter0,rpt));
    t_map = Image(getSnapshotImage(figReporter0,rpt));

    t_map.Width = '4.7in';
    t_map.Height = '4.1in';

    add(rpt,t_map);

    saveas(gcf, ([figureDir (sprintf('/%s_maps.png',  mapsList{thisFig}))]));
    
     


end

% tab = Table({t_map{1},t_map{2},t_map{3},t_map{4};mapsList{1},mapsList{2},mapsList{3},mapsList{4}})
% add(rpt,tab);



%%

chapter = Chapter();

ct = 1;


for r = 1 : length(path2roi)
    for hemi = 1 :length(hemispheres)
        
        
        
        roi = [];
        
        ind  = read_label(['sub-' subject],sprintf ('%s.%s%s',hemispheres{hemi},path2roi{r}));
        roi  = [roi; ind(:,1) + 1];
        myroi = zeros(size(map_file.vexpl.(hemispheres{hemi})));
        myroi(roi) = 1;
        
        %         thr  = double(map_file.vexpl.(hemispheres{hemi})>0.15) & double(map_file.eccen.(hemispheres{hemi})<10) & myroi & double(map_file.sigma.(hemispheres{hemi})>0.25);
        thr  = double(map_file.vexpl.(hemispheres{hemi})>0.2) & double(map_file.eccen.(hemispheres{hemi})<12) & double(map_file.eccen.(hemispheres{hemi})>0) & double(map_file.sigma.(hemispheres{hemi})>0) & myroi;
        
        
        fig2=figure(2);
        sgtitle(sprintf('Var thr > %.1f, Ecc < %.0f',0.2,12))
        subplot(2,2,ct)
        prfsize = map_file.sigma.(hemispheres{hemi})(thr);
        ecc = map_file.eccen.(hemispheres{hemi})(thr);
        


        scatter(ecc, prfsize,'.','k')
        xlabel('Eccentricity');
        ylabel('pRF size');
        R = corrcoef(ecc,prfsize);
        % add line of best fit
        coeff = polyfit(ecc, prfsize, 1);
        yFit = polyval(coeff , 0.2:0.1:12);
        hold on;
        linearfit= plot(0.2:0.1:12, yFit, 'b-', 'LineWidth', 2.5);
        
        grid on;
        title(sprintf('size vs eccen, %s %s',path2roi{r},hemispheres{hemi}),'Interpreter','None')
        xlim([0 12])
        ylim([0 6])
        legend(linearfit,['y = ' num2str(coeff(1)) '*x + ' num2str(coeff(2))]);
        ct = ct + 1;
        %         title(sprintf('%s %s',path2roi{r},hemispheres{hemi}))
        set(gca,'Fontsize',10)
        set(gcf,'Position',[179   313   486   496])
    end
    
end

saveas(gcf, ([figureDir (sprintf('/%s_size_vs_eccen.png',  mapsList{thisFig}))]));

figReporter0 = Figure(fig2);
chapter.Title = 'pRFplots';
add(rpt,chapter);
figReporter0.Snapshot.Height = [];
figReporter0.Snapshot.Width = [];
add(rpt,figReporter0);

rptview(rpt)
%% relevant sub functions
    function plot_mesh(faces, vertices, map, cmap)
        t = trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), map, 'FaceColor', 'flat');
        t.LineStyle = 'none';
        axis equal; hold on;
        %     cmap = [[1 1 1]*.7; cmap];
        colormap(gcf,cmap);
    end

    function set_view(gcf)
        axis off; set(gcf, 'color','white','InvertHardCopy', 'off');
        view(-45,0);
        material dull;
        h=light; lightangle(h,  45, 45); lighting gouraud;
        h=light; lightangle(h, -45, 45); lighting gouraud;
        h=light; lightangle(h, -45, -90); lighting gouraud;
        set(gcf, 'Position', [150 100 750 625]);
        axis tight
    end

end

%% ******************************
% ******** SUBROUTINES **********
% *******************************

function latestDir = find_latest_dir(projectDir)


d = dir(sprintf('%s/derivatives/',projectDir));
d = d([d(:).isdir]==1);
[~,id] = sort([d.datenum]);
d = d(id);
latestDir = d(end).name;

end