%% PALA_SilicoFlow_fig.m : Post Processing - errors and pairings algorithms for IN SILICO FLOW
% load results from localization and pairings algorithms and create displaying figures.
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
%% Select IQ and Media folder
fprintf('Running PALA_SilicoFlow_fig.m - Images generation\n')
workingdir = [PALA_data_folder '\PALA_data_InSilicoFlow'];cd(workingdir)
filename = 'PALA_InSilicoFlow';

myfilepath = [workingdir filesep filename];
myfilepath_data = [workingdir filesep 'Results' filesep filename];
myfilepath_fig = [workingdir filesep 'Results' filesep 'img' filesep];

listVar = {'P','PData','Trans','Media','UF','Resource','Receive','filetitle'};
load([myfilepath '_sequence.mat'],'-mat',listVar{:});clear listVar

MatOutTarget = load([workingdir filesep filename '_v3_config.mat'],'MatOut');MatOutTarget = MatOutTarget.MatOut;
load([myfilepath_data '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')
load([myfilepath_data '_MatOut_multi_' num2str(abs(30)) 'dB'  '.mat'],'-mat')

if exist(myfilepath_fig)~=7;mkdir(myfilepath_fig);end
[listAlgoName,ListColor,ListMarker,ListShortName] = PALA_GetFormat; labelSize = 15;

%% FIG 5: MATOUT target %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig39 = figure(39);clf;fig39.Position = [20 150 size(MatOutTarget,2)/2 size(MatOutTarget,1)/2];fig39.InvertHardcopy='off';

tight_subplot(1,1,0,0,0.0);
im = imagesc(imgaussfilt((MatOutTarget).^(1/2),1));
im.CData = min(im.CData,10); % caxis
axis image, colormap hot, axis off

posScale(1) = 60; posScale(2) = size(MatOutTarget,1)-70;
hold on;caxis([0 10])
im.CData(posScale(2)+[1:10],posScale(1)+[1:100])=max(caxis);
tscale = text(posScale(1)+50, posScale(2),'1 mm','color','w','fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center');

figname = 'fig5_MatOutTarget';
print([myfilepath_fig figname],'-dpng','-r300')
WriteTif(im.CData,hot(128),[myfilepath_fig figname '.tif'],'overwrite',1)

%%
fig = figure(38);clf;fig.Units = 'centimeters';fig.Position = [7 2 1 8];
ipow = 1/2;caxis([0 10]);colormap(hot(128))
clb = colorbar('Position',[.01 .05 .1 .9],'AxisLocation','in');
val = 0:10:(10^(1/ipow));clb.Ticks = val.^ipow;
clb.TickLabels = val;clb.TickLabels(2:2:end,:)=' ';clb.TickLength =0.02;clb.Label.String='Counts';
clb.Label.Position = [1 2];clb.FontSize = labelSize-4;axis off
figname = 'fig5_MatOutTarget_clb';
print([myfilepath_fig figname],'-dpng','-r300')
print(fig,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% FIG 5: MATOUTS zoom horseshoe %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig40 = figure(40);clf;fig40.Position = [169 511 688 426];fig40.InvertHardcopy='off';
a = tight_subplot(2,3,[0.02 0.005],0,0.0);

list2Disp = [1 2 3 6 7 8];
xcrop = [230:360]; ycrop = [210:330];
figname = 'fig5_MatOut_horseshoe_30dB';
for ii=1:numel(list2Disp)
    ialgo = list2Disp(ii);axes(a(ii))
    if ii==1
        im = imagesc(imgaussfilt((MatOutTarget(ycrop,xcrop)).^(1/2),.001));
        tname(ialgo) = text(3, 3,'Target MatOut','color','w','fontsize',labelSize-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    else
        im = imagesc(imgaussfilt((MatOut{ialgo-1}(ycrop,xcrop)).^(1/2),.001));
        tname(ialgo) = text(3, 3,listAlgoName{ialgo-1},'color','w','fontsize',labelSize-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    end
    im.CData = min(im.CData,10);caxis([0 10]) % caxis
    axis image, colormap hot;axis off
    if ii==1
        posScale(1) = 8; posScale(2) = size(im.CData,1)-5;
        hold on
        im.CData(posScale(2)+[0:2],posScale(1)+[1:50])=max(caxis);
        tscale = text(posScale(1)+25, posScale(2),'500 \mum','color','w','fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center','Interpreter','tex');
    end
    if ii==1,WriteTif(im.CData,hot(128),[myfilepath_fig figname 'Target' '.tif'],'overwrite',1)
    else,WriteTif(im.CData,hot(128),[myfilepath_fig figname ListShortName{ialgo-1} '.tif'],'overwrite',1)
    end
end
print([myfilepath_fig figname],'-dpng','-r300') % Save fig

%% FIG 5: MATOUTS zoom 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig42 = figure(42);clf;fig42.Position = [169 257 1133 451];fig42.InvertHardcopy='off';
a = tight_subplot(2,3,[0.01 0.01],0,0.0);
MatZoom.x = 250:650;MatZoom.z = 350:590;

figname = 'fig5_MatOut_zoom_30dB';
for ii=1:numel(list2Disp)
    ialgo = list2Disp(ii);axes(a(ii))
    if ii ==1
        im = imagesc(imgaussfilt((MatOutTarget(MatZoom.z,MatZoom.x)).^(1/2),.5));
        tname(ialgo) = text(5, 5,'Target MatOut','color','w','fontsize',labelSize+2-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    else
        im = imagesc(imgaussfilt((MatOut{ialgo-1}(MatZoom.z,MatZoom.x)).^(1/2),.001));
        tname(ialgo) = text(5, 5,listAlgoName{ialgo-1},'color','w','fontsize',labelSize+2-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    end
    im.CData = min(im.CData,10); % caxis
    axis image, colormap hot,axis off
    if ii==1
        posScale(1) = 10; posScale(2) = size(im.CData,1)-5;
        im.CData(posScale(2)+[0:3],posScale(1)+[1:50])=max(caxis);
        hold on
        tscale = text(posScale(1)+25, posScale(2),'500 \mum','color','w','fontsize',labelSize+2-2,'VerticalAlignment','bot','HorizontalAlignment','center','Interpreter','tex');
    end
    if ii==1,WriteTif(im.CData,hot(128),[myfilepath_fig figname 'Target' '.tif'],'overwrite',1)
    else,WriteTif(im.CData,hot(128),[myfilepath_fig figname ListShortName{ialgo-1} '.tif'],'overwrite',1)
    end
end
print([myfilepath_fig figname],'-dpng','-r300') % Save fig

%% FIG 3: Error lat ax flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([myfilepath_data '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')
fig43 = figure(20);clf;fig43.Position = [223 414 1250 340];

aa0 = tight_subplot(2,7,[.09 .015],[.07 .01],[.01 .01]);
xedge = linspace(-.55,.55,200);

labelSize = 16;
for ii = 1:numel(listAlgo)
    % for each algo display histograms
    ialgo=ii;axes(aa0(ii));
    axhist_x = histogram((ErrList{ialgo}(:,3)),xedge,'edgecolor','none','Normalization','Probability','FaceColor',[1 1 1]*0.2);
    axhist_x.FaceAlpha = 1;axhist_x.FaceColor = ListColor(ii,:);
    ErrMean_x(ialgo) = mean(ErrList{ialgo}(:,3));
    ErrStd_x(ialgo) = std(ErrList{ialgo}(:,3));
    hold on;xlim([-1 1]*xedge(end)),ylim([0 .06]);grid on

    text(.02,aa0(ii).YLim(2)-.009,['\sigma=' num2str(round(ErrStd_x(ialgo),2),'%.2f') ],'fontsize',labelSize,'HorizontalAlignment','left','VerticalAlignment','top')

    ax = gca;ax.Color = 'none';
    ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';
    ax.YTickLabel = [];ax.XTickLabel = {'-\lambda/2','','0','','\lambda/2'};ax.XTick = [-.5,-.25,0,.25,.5];
    ax.FontSize = labelSize-4;
    ax.TickDir ='in';ax.TickLength = [.04 .2];ax.Box='off';
    ax.GridAlpha = .5; ax.GridLineStyle = '--';
%     title(listAlgoName{ialgo},'color',cfg.title,'fontsize',labelSize,'interpreter','none')
    if ii==1;text(-.55,mean(ylim),'Lateral error','FontSize',labelSize,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90);end
end

for ii= 1:numel(listAlgo)
    ialgo=ii;axes(aa0(ii+7));
    axhist_z = histogram(ErrList{ialgo}(:,2),xedge,'edgecolor','none','Normalization','Probability','FaceColor',[1 1 1]*0.2);
    axhist_z.FaceAlpha = axhist_x.FaceAlpha;
    ErrMean_z(ialgo) = mean(ErrList{ialgo}(:,2));
    ErrStd_z(ialgo) = std(ErrList{ialgo}(:,2));
    axhist_z.FaceColor = ListColor(ii,:);
    hold on;xlim([-1 1]*xedge(end));ylim([0 .03]);grid on
    text(.02,aa0(ii+7).YLim(2)-.001,['\sigma=' num2str(round(ErrStd_z(ialgo),2),'%.2f') ],'fontsize',labelSize,'HorizontalAlignment','left','VerticalAlignment','top')

    set(gca, 'yticklabel', []);set(gca, 'xticklabel', {'-\lambda','0','+\lambda'});
    az = gca; az.Color = ax.Color;
    az.XAxisLocation = 'origin';az.YAxisLocation = 'origin';
    az.YTickLabel = [];az.XTickLabel = ax.XTickLabel;az.XTick = ax.XTick;
    az.FontSize = ax.FontSize;az.TickDir = ax.TickDir;az.TickLength = ax.TickLength;az.Box= ax.Box;
    az.GridAlpha = ax.GridAlpha; az.GridLineStyle = ax.GridLineStyle;
    if ii ==1,text(-.55,mean(ylim),'Axial error','FontSize',labelSize,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90);end
end
figname = 'fig3_distrib_30dB';
% print(fig43,[myfilepath_fig figname],'-dpng','-r300');% Save fig
% print(fig43,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%%
% ErrMean_x_interp = mean([ErrList{3}(:,3);ErrList{4}(:,3);ErrList{5}(:,3)]);
% ErrStd_x_interp = std([ErrList{3}(:,3);ErrList{4}(:,3);ErrList{5}(:,3)]);
% ErrMean_z_interp = mean([ErrList{3}(:,2);ErrList{4}(:,2);ErrList{5}(:,2)]);
% ErrStd_z_interp = std([ErrList{3}(:,2);ErrList{4}(:,2);ErrList{5}(:,2)]);

%% FIG 3: TRUE POSITIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([myfilepath_data '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')
TP = squeeze(Stat_class(3,:,:));FN = squeeze(Stat_class(4,:,:));FP = squeeze(Stat_class(5,:,:));

J_p = TP./(FP+TP)*100;      % precision
J_r = TP./(FN+TP)*100;      % sensitivity (recall)
J_ac = TP./(FP+FN+TP)*100;  % Jaccard
J_p_mean = mean(J_p,2);J_r_mean = mean(J_r,2);J_ac_mean = mean(J_ac,2);
TP = sum(TP,2);FN = sum(FN,2);FP = sum(FP,2);

%%
fig50=figure(50);clf;hold on
fig50.Position = [880 600 708 250];
b1(1) = barh(1:7,TP*1e-3,1,'stack','DisplayName','True Positive');
b1(2:3) = barh(1:7,[-FN,-FP]*1e-3,1,'stack');

b1(1).FaceColor = 'flat';b1(1).CData = ListColor;b1(1).FaceAlpha =1;
b1(2).FaceColor = 'flat';b1(2).CData = ListColor;b1(2).FaceAlpha =.5;b1(2).DisplayName = 'False Negative';
b1(3).FaceColor = 'flat';b1(3).CData = ListColor;b1(3).FaceAlpha =.2;b1(3).DisplayName = 'False Positive';
legend('box','off','Location','northwest')

for ii=1:numel(listAlgo)
   text(TP(ii)*1e-3+2e1,ii,listAlgoName(ii))
   text(TP(ii)*1e-3/2,ii, [num2str(round(TP(ii)/1000)) 'k'],'HorizontalAlignment','center')
   text(-FN(ii)*1e-3/2,ii, [num2str(round(FN(ii)/1000)) 'k'],'HorizontalAlignment','center')
   text(-FP(ii)*1e-3/2-FN(ii)*1e-3,ii, [num2str(round(FP(ii)/1000)) 'k'],'HorizontalAlignment','center')
end
xlim([-7*1e2 5*1e2]);ylim([.5 7.5])
a=gca;a.YTick=[];a.YAxisLocation = 'origin';a.XAxis.Exponent = 0;xtickformat('%.0fk')
a.XAxis.TickDirection = 'out';
a.Position = [.05 .12 .9 .87];

text(a.XLim(1),.5,'False','HorizontalAlignment','right','VerticalAlignment','middle')
text(a.XLim(2),.5,'True','HorizontalAlignment','left','VerticalAlignment','middle')

NbDetectableBulles = mean(TP+FN);
plot([1 1]*NbDetectableBulles/1000,[0 7.5],'--','LineWidth',1,'color',[0 0 0 .2],'HandleVisibility','off');grid on
text(NbDetectableBulles/1000,.2,'Max','HorizontalAlignment','center')

%% Save fig
figname = 'fig3_TP_FP_FN_30dB';
print(fig50,[myfilepath_fig figname],'-dpng','-r300')
print(fig50,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% FIG 3: JACCARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig52 = figure(52);clf;fig52.Position = [223 590 563 244];
barcolor = [1 1 1]*.5;stsuba = tight_subplot(1,3,[0 .02],[.15 .1],[.01 .05]);

axes(stsuba(1));
b_p = barh(J_p_mean,'facecolor',barcolor,'BarWidth',1);
text(b_p.XData,b_p.YData+20,string(round(b_p.YData,1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',labelSize)
b_p.FaceColor = 'flat';
b_p.CData = ListColor;
text(J_p_mean-1,1:length(J_p_mean),num2str(round(J_p_mean,1)),'VerticalAlignment','middle','HorizontalAlignment','right');
hold on,errorbar(J_p_mean,1:7,std(J_p,0,2),'horizontal','Color','k','LineStyle','none','CapSize',10);

ylim([0.5 7.5]);xlim([0 90])
ap = gca;ap.XAxisLocation = 'origin'; ap.YAxisLocation = 'origin';
ap.YTickLabel = {};ap.FontSize = labelSize-2;ap.TickDir = 'in';ap.Box = 'off';
ap.XAxis.FontSize = labelSize -4;grid on;xticks([0 30 60])
text(5,7.5,'Precision','FontSize',labelSize,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold')

axes(stsuba(2));suba = gca;
b_r = barh(J_r_mean,'facecolor',barcolor,'BarWidth',1);
text(b_r.XData,b_r.YData+20,string(round(b_r.YData,1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',labelSize)
hold on,errorbar(J_r_mean,1:7,std(J_r,0,2),'horizontal','Color','k','LineStyle','none','CapSize',10);
b_r.FaceColor = 'flat';
b_r.CData = ListColor;
text(J_r_mean-1,1:length(J_r_mean),num2str(round(J_r_mean,1)),'VerticalAlignment','middle','HorizontalAlignment','right');

ylim([0.5 7.5]);xlim([0 80])
ar = gca;ar.XAxisLocation = 'origin';ar.YAxisLocation = 'origin';
ar.YTickLabel = {};ar.FontSize = ap.FontSize;ar.TickDir = ap.TickDir;
ar.XAxis.FontSize = ap.XAxis.FontSize;ar.Box= 'off';grid on;xticks([0 30 60])
text(5,7.5,'Sensitivity','FontSize',labelSize,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold')

axes(stsuba(3));suba = gca;
b_ac = barh(J_ac_mean,'facecolor',barcolor,'BarWidth',1);
text(b_ac.XData,b_ac.YData+20,string(round(b_ac.YData,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',labelSize)
hold on,errorbar(J_ac_mean,1:7,std(J_ac,0,2),'horizontal','Color','k','LineStyle','none','CapSize',10);
b_ac.FaceColor = 'flat';
b_ac.CData = ListColor;
text(J_ac_mean-1,1:length(J_ac_mean),num2str(round(J_ac_mean,1)),'VerticalAlignment','middle','HorizontalAlignment','right');

ylim([0.5 7.5]);xlim([0 100])
aa = gca;aa.XAxisLocation = 'origin';aa.YAxisLocation = 'origin';
aa.FontSize = ap.FontSize;aa.TickDir = ap.TickDir;
aa.XAxis.FontSize = ap.XAxis.FontSize;aa.YTickLabel = {};
aa.Box = 'off';grid on;xticks([0 30 60])

text(b_ac.YData+5,b_ac.XData,listAlgoName(1:7),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize-2)
text(5,7.5,'Jaccard Index','FontSize',labelSize,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold')

%% Save fig
figname = 'fig3_Jacc_30dB';
print(fig52,[myfilepath_fig figname],'-dpng','-r300')
print(fig52,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% %%%%%%% FIG 2 : BRANCHES MATOUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MM = MatOutTarget(780:800,120:800);DistNul = [];Gap = [];ListGap = [];
Threshold_p = .03;
Threshold_val = Threshold_p.*max(MM(:,1));
for ialgo = 1:8
    if ialgo<8, MM = MatOut{ialgo}(780:800,220:800);else, MM = MatOutTarget(780:800,220:800);end
    lz = 1:2:size(MM,2);
    for ind=1:numel(lz)
        Lseq = MM(:,lz(ind));Lseq(Lseq<0)=0;
        Lseq = Lseq(find(Lseq>Threshold_val,1,'first'):find(Lseq>Threshold_val,1,'last'));
        DistNul(ind,ialgo) = sum(Lseq<Threshold_val);
        if ialgo==8
            im = imregionalmax(Lseq);
            posMax = find(im); Kepmax = [0;diff(posMax)]~=1; posMax=posMax(Kepmax);
            val = Lseq(im);
            [~,b2] = sort(val(Kepmax),'descend');
            try Gap(ind) = abs(posMax(b2(1))-posMax(b2(2)));catch,Gap(ind)=0;end
        end
    end
    [mingap,~] = find(DistNul(:,ialgo)==0,1,'last');
    ListGap(ialgo) = mingap+1;
end
SR2conv = 1./ULM.res; % convert in \lambda
clear mingap ind ialgo MM im posMax Kepmax b2 val
%%
DistNul_wv = DistNul*SR2conv; clear DistNul
fig26 = figure(26);fig26.Position = [500 200 1052 493];clf; clear plt_dist
ax1 = axes('Position',[.05 .35 .92 .6]);hold on

Gap = DistNul_wv(:,end); % take the measured gap on the Traget MatOut in [nbSRpixel]
pp = polyfit(lz(:),Gap(:),1);
Gap_linear = pp(2) + lz*pp(1);
ListGap_pix = Gap_linear(ListGap);

for ialgo=1:7
    plt_dist(ialgo) = plot(Gap_linear(:),DistNul_wv(:,ialgo),'Marker',ListMarker{ialgo},'MarkerSize',7,'DisplayName',listAlgoName{ialgo},'color',ListColor(ialgo,:),'linewidth',1.5);
    plt_dist(ialgo).MarkerFaceColor = plt_dist(ialgo).Color;
    plt_dist(ialgo).MarkerIndices = plt_dist(ialgo).MarkerIndices((ialgo*3):20:end);
    plot(ListGap_pix(ialgo),0,'+','color',ListColor(ialgo,:),'MarkerSize',10)
    text(ListGap_pix(ialgo),0,num2str(round(ListGap_pix(ialgo),2)),'VerticalAlignment','top','HorizontalAlignment','center')
end

lgd = legend([plt_dist],'Location','northwest','Interpreter','tex');
xlim([3 9]*SR2conv);ylim([0 12]*SR2conv);grid on
ylabel('Measured canal to canal distance [\lambda]');xlabel('Simulated canal to canal distance')
title(['Size of the measured gap on Matout (threshold ' num2str(Threshold_val) ', ' num2str(Threshold_p*100) '%)'])

figname = 'fig2_gap_30dB';
for ialgo=0:7
    ax(ialgo+1) = axes('Position',[.05 + (ialgo)/7*.8,.01,.8/8,.2]);
    if ialgo>0,MM = MatOut{ialgo}(780:800,220:800);
    else,MM = MatOutTarget(780:800,220:800);end
    im=imagesc(MM.^(1/2));axis off, colormap hot;aa = gca;
    if ialgo>0,title(listAlgoName{ialgo},'Interpreter','none');savName = ListShortName{ialgo};else title('Target','Interpreter','none');savName='TG';end
    WriteTif(im.CData,aa.Colormap,[myfilepath_fig figname savName '.tif'],'caxis',caxis,'overwrite',1)
end

%% Save fig
figname = 'fig2_gap_30dB';
print(fig26,[myfilepath_fig figname],'-dpng','-r300')
print(fig26,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')
clear MM im DistNul_wv Gap

%% Fig 2 : RMSE(dB) silico %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([myfilepathdata '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')
val_Err = {};RMSE_m=[];
for ii=1:numel(ClutterList)
    temp = load([myfilepath_data '_Stats_multi' num2str(abs(ClutterList(ii))) 'dB' '.mat'],'ErrList');
    for ialgo=1:numel(listAlgo)
        val_Err{ialgo,ii} = temp.ErrList{ialgo}(:,1);
        if ClutterList(ii)==-30;RMSE_m(ialgo) = mean(val_Err{ialgo,ii});end
    end
end
clear temp

%% build figure
fig60 = figure(60);clf;fig60.Position = [190 470 1150 240];
axes('Position',[.04 .05 .96 .95])
vv = PALA_boxplot_multi(listAlgoName,val_Err,15,ListColor',10,0);ylim([0 .56])

for ii=1:numel(ClutterList)
    tt=text(.55+(ii-1)/(numel(ClutterList)-1)*.8,.01,[num2str(abs(ClutterList(ii)))],'VerticalAlignment','bot','HorizontalAlignment','left','FontSize',9);
end
tt.String = [tt.String ' dB'];grid on
a = gca;a.FontSize = 11;a.YTick = [0:.1:.50];
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;
%% Save fig
figname = 'fig2_RMSEsilico';
print(fig60,[myfilepath_fig figname],'-dpng','-r300')
print(fig60,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% 2D plot RMSE SNR
fig61 = figure(61);clf;fig61.Position = [706 601 520 341];hold on
for ialgo=1:Nalgo
    plot(abs(ClutterList),squeeze(vv(ialgo,:,1)),'.-','DisplayName',listAlgoName{ialgo},...
        'Color',ListColor(ialgo,:),'Marker',ListMarker{ialgo},'LineWidth',1.5,'MarkerSize',10,'MarkerFaceColor',ListColor(ialgo,:),'MarkerEdgeColor','none')
end
ylabel('RMSE');ylim([0.06 .35]);xlabel('SNR [dB]');grid on
ca = gca;ca.FontSize = 11;ca.XDir = 'reverse';
legend;ca.Legend.Location = 'south';
%% Save fig
figname = 'fig2_RMSEsilicoSNR';
print(fig61,[myfilepath_fig figname],'-dpng','-r300')
print(fig61,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% Save data for Global Score
Radar.RMSE          = RMSE_m;
Radar.Jaccard       = squeeze(mean(J_ac,2));
Radar.precision     = squeeze(mean(J_p,2));
Radar.Gap           = ListGap_pix;
Radar.listName      = listAlgoName;
save([workingdir filesep filename '_scores'],'Radar','Nalgo')

fprintf('PALA_SilicoFlow_fig.m done.\n');

