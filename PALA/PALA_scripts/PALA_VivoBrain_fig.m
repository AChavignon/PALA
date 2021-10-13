%% PALA_VivoBrain_fig.m : Post Processing - displaying and anylisng results for IN VIVO BRAIN
% For each algorithm, bubbles are detected, localised and tracks.
% This script compare results of differents algorithms building matout, and comparing various behaviours.
%
% Created by Arthur Chavignon 25/02/2020
%
% DATE 2020.12.17 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

run('PALA_SetUpPaths.m')
%% Selected data file and saving folders
fprintf('Running PALA_VivoBrain_fig.m - Images generation\n')
workingdir = [PALA_data_folder '\PALA_data_InVivoRatBrain'];cd(workingdir)
filename = 'PALA_InVivoRatBrain_';

SavingFolder = 'Results';
myfilepath = [workingdir filesep filename];
myfilepath_data = [workingdir filesep SavingFolder filesep filename];
myfilepath_fig = [workingdir filesep SavingFolder filesep];
if exist(myfilepath_fig)~=7;mkdir(myfilepath_fig);end

load([myfilepath_data 'MatOut_multi.mat'])
load([myfilepath_data 'MatOut_multi_nointerp.mat'],'ProcessingTime','ULM','listAlgo','MatOutNoInterp')
[listAlgoName,ListColor,ListMarker,ListShortName] = PALA_GetFormat;labelSize = 15;Nalgo = numel(listAlgo);clear MatOut_vel
%% %% FIG 6 : MATOUT RAIDAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig38 = figure(38);clf;fig38.Position = [5 160 1900 830];
countpos = [-20 -20];

ialgo = 7; % select radial symmetry algorithm
ii=imagesc(imgaussfilt(MatOut{ialgo}(17:765,23:1166).^(1/3),.01));

fig38.Position(4:-1:3)=size(ii.CData);fig38.InvertHardcopy = 'off';
aa = gca;
aa.Position = [0 0 1 1];
axis image, colormap(hot(128)),axis off
text(12,10,listAlgoName{ialgo},'color','w','fontsize',labelSize,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');

posScale(1) = 60; posScale(2) = size(MatOut{ialgo},1)-70;hold on
pscale = plot(posScale(1)+[1 100],posScale(2).*[1 1],'w-','linewidth',3);
tscale = text(mean(pscale.XData), pscale.YData(1),'1mm','color','w','fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center');
countpos = size(ii.CData)+countpos;
ii.CData(posScale(2)+[0:3],posScale(1)+[1:100])=max(caxis);
ii.CData = min(ii.CData,10);caxis([0 10])

tcount = text(countpos(2), countpos(1),[num2str(round(NbrOfLoc(ialgo)/1e6,2)),'M'],'color','w','fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','right');
figname = 'fig6_matout_radial_full';
print([myfilepath_fig figname ],'-dpng','-r300')
WriteTif(ii.CData,aa.Colormap,[myfilepath_fig figname '.tif'],'caxis',caxis,'overwrite',1)

%%
fig = figure(38);clf;fig.Units = 'centimeters';fig.Position = [7 2 8 1];
ipow = 1/3;caxis([0 10]);colormap(hot(128))
clb = colorbar('Location','south','Position',[.04 .01 .9 .3]);
val = [0:10:100 200:100:1000];
clb.Ticks = val.^ipow;
clb.TickLabels = val;clb.TickLabels([3:5 7:10 13:14 16:19],:)=' ';
clb.AxisLocation = 'in';clb.Label.String='Counts';
clb.TickLength = .02;
clb.Label.Position = [.9 1.5];clb.FontSize = labelSize-4;axis off
figname = ['fig6_matout_radial_full_clb'];drawnow
print([myfilepath_fig figname ],'-dpng','-r300')
print(fig,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% %% FIG 6 : MATOUTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig40 = figure(40);clf;fig40.Position = [5 160 1290 860];fig38.InvertHardcopy = 'off';
MatZoom.x = 476:627;MatZoom.z = 651:882;
a = tight_subplot(2,2,[0.01 0.005],0.00,0.00);
countpos = [-20 -20];
figname = 'fig6_matout_radial_full';
listAlgoMatOut = [2 3 5 7];
for ii=1:numel(listAlgoMatOut)
    ialgo=listAlgoMatOut(ii);
    axes(a(ii))
    im=imagesc((MatOut{ialgo}(MatZoom.x,MatZoom.z)).^(1/3));
    axis image, colormap hot, axis off
    im.CData = min(im.CData,10);caxis([0 10]) % caxis 10

    if ii==1
        posScale(1) = 10;
        posScale(2) = round(a(ii).YLim(2)-10);
        hold on
        pscale = plot(posScale(1)+[1 50],posScale(2).*[1 1],'w-','linewidth',3);
        tscale = text(mean(pscale.XData), pscale.YData(1),'500\mum','color','w','fontsize',labelSize+4,'VerticalAlignment','bot','HorizontalAlignment','center','Interpreter','tex');
        countpos = [a(ii).YLim(2) a(ii).XLim(2)]+[-5,-10];
        im.CData(posScale(2)+[0:2],posScale(1)+[1:50])=max(caxis);
    end
    text(2,2,listAlgoName{ialgo},'color','w','fontsize',labelSize+4,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    tcount = text(countpos(2), countpos(1),[num2str(round(NbrOfLoc(ialgo)/1e6,2)),'M'],'color','w','fontsize',labelSize+2,'VerticalAlignment','bot','HorizontalAlignment','right');
    WriteTif(im.CData,hot,[myfilepath_fig figname num2str(ii) '.tif'],'caxis',caxis,'overwrite',1)
end
axes(a(end)),axis off,linkaxes(a)

figname = 'fig6_matouts_tight';
print([myfilepath_fig figname],'-dpng','-r300') % Save fig

%% %% FIG 6 : GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Gridding Index with MatOut WITHOUT interpolation
[GridIndex,F_0] =  PALA_ComputeGriddingIndex(MatOutNoInterp);

fig41 = figure(41);clf;fig41.Position = [530 410 310 220];
bh = barh(GridIndex,'facecolor',[1 1 1]*0.5,'BarWidth',1);
bh.FaceColor = 'flat';
bh.CData = ListColor;

text(bh.YData-.1,bh.XData,string(round(GridIndex,1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',labelSize-2)
text(bh.YData+.5,bh.XData,listAlgoName,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize)

ae = fig41.Children;
ylim([0 8] + [1 -1]*.5);xlim([0 23]);hold on

ae.Position = [.02 .13 .8 .75];
ae.XAxisLocation = 'origin'; ae.YAxisLocation = 'origin';
ae.FontSize = labelSize-2;ae.XAxis.FontSize = labelSize;
ae.TickDir = 'both';ae.TickLength = [.005 .0];
ae.Box= 'off';ae.GridAlpha = .4;ae.GridLineStyle = '--';
ae.XLabel.String = '[AU]';ae.XLabel.Position(1:2)= [23 0.3];
ae.XAxis.FontSize = labelSize-2;
set(gca,'ytick',[])
set(get(ae,'children'),'clipping','off')% turn off clippings
title('Gridding Index','FontSize',labelSize)

figname = 'fig6_GriddingIndex';
print(fig41,[myfilepath_fig figname],'-r300','-dpng') % Save fig
print(fig41,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% %% FIG 6 : Saturation curves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure(78);fig.Position = [841 526 673 359];
lz = [1:size(MatOutSat,1)]*800/1000;
SaturationEnd = MatOutSat(end,:)/numel(MatOut{1});

plt = plot(lz,MatOutSat/numel(MatOut{1})*100,'LineWidth',1.2);
marSifht = [0 0 10 20 30 13 26];
for ii=1:numel(plt)
    plt(ii).Color = ListColor(ii,:);
    plt(ii).Marker = ListMarker{ii};
    plt(ii).MarkerFaceColor = ListColor(ii,:);
    plt(ii).MarkerIndices = [marSifht(ii)+10:50:size(MatOutSat,1)];
    plt(ii).MarkerSize = 10;
end
legend(listAlgoName,'Location','best','FontSize',13)
xlabel('Number of frames')
ylim([0 50]);xticks([0:25:200]);xtickformat('%,.0fk');

grid on;a = gca;a.Position = [.07 .12 .9 .8];a.FontSize=labelSize;
a.XTickLabel{2}='';a.XTickLabel{4}='';a.XTickLabel{6}='';a.XTickLabel{8}='';
title('Saturation (%)')

%%
figname = 'fig6_saturation';
print(fig,[myfilepath_fig figname],'-r300','-dpng')
print(fig,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% %% FIG 2 : Processing Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig45 = figure(45);clf;fig45.Position = [530 400 390 210];
ae = axes('Position',[.01 .11 .92 .75]);
TotalProTime = mean(ProcessingTime,1)/60;TotalProTime = TotalProTime./min(TotalProTime);

b_ap = barh(TotalProTime,'facecolor',[1 1 1]*0.5,'BarWidth',1);
hold on,b_ap_err = errorbar(TotalProTime,1:Nalgo,std(ProcessingTime,0,1)/60/min(TotalProTime),'horizontal','Color','k','LineStyle','none','CapSize',10);
b_ap.YData(end-1)=7.5;b_ap_err.XData(end-1)=b_ap.YData(end-1);
b_ap.FaceColor = 'flat';b_ap.CData = ListColor;

text(b_ap.YData+.1,b_ap.XData,string(round(TotalProTime,2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize-2)
grid on
ylim([0 8] + [1 -1]*.5)

ae.XAxisLocation = 'origin'; ae.YAxisLocation = 'origin';
ae.FontSize = labelSize-2;ae.XAxis.FontSize = labelSize;
ae.TickDir = 'both';ae.TickLength = [.005 .0];ae.Box= 'off';
ae.XLabel.Position(1:2)= [23 0.3];
ae.XAxis.FontSize = labelSize-4;
set(gca,'ytick',[]);ae.XTick = [1:7];
title('Processing time','FontSize',labelSize)

figname = 'fig2_processingTime';
print(fig45,[myfilepath_fig figname],'-r300','-dpng')
print(fig45,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% %% Intenstity MatOutnointerp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nloc = cellfun(@(x) sum(x(:)),MatOutNoInterp);
TotalProTime = sum(ProcessingTime,1);

%% RADAR DATA
Radar_vivo.Saturation = SaturationEnd;
Radar_vivo.GridIndex = GridIndex;
Radar_vivo.Time = TotalProTime;
Radar.listName = listAlgoName;
save([workingdir filesep filename 'scores'],'Radar_vivo')

fprintf('PALA_VivoBrain_fig.m done.\n');
