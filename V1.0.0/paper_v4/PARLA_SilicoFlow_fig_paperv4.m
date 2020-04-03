%% ULM PARLA: Post Processing - errors and pairings algorithms for IN SILICO FLOW
% load results from localization and pairings algorithms and create displaying figures.
%
% Created by Arthur Chavignon 25/02/2020
%
% DATE 2020.03.30 - VERSION 1.0.0
% AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
% Code Available under Creative Commons Non-Commercial 4.0
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018 

clear all;close('all')
PARLA_addons_folder = 'F:\ArthurC\CHAPRO_local\_SIMULATIONS\FightClub_PostPro\addons'; % location of the addons folder
SimppleTracker_folder = 'F:\ArthurC\CHAPRO_local\_SIMULATIONS\FightClub_PostPro\SimpleTracker'; % location of the SimpleTracker
PARLA_data_folder = 'F:\ArthurC\Data\Fight_Club\Data'; % path of data
addpath(genpath(PARLA_addons_folder))
addpath(genpath(SimppleTracker_folder))
%% Select IQ and Media folder
workingdir = [PARLA_data_folder '\PARLA_data_InSilicoFlow'];
filename = '2D_flow_1203-180447';
cd(workingdir)

myfilepath = [workingdir filesep filename];
myfilepathdata = [workingdir filesep 'config2' filesep filename];

listVar = {'P','PData','Trans','Media','UF','Resource','Receive','filetitle'};
load([myfilepath '_sequence.mat'],'-mat',listVar{:})

suffixe = '2';

MatOutTarget = load([workingdir filesep '2Dpaper3_config.mat'],'MatOut');MatOutTarget = MatOutTarget.MatOut;
% load([myfilepath '_Tracks_multi.mat'],'-mat')
load([myfilepathdata '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')
load([myfilepathdata '_MatOut_multi_' num2str(abs(30)) 'dB'  '.mat'],'-mat')

[listAlgoName,ListColor,ListMarker,ShortName] = GetFormatFightClub;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE MATOUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelSize = 15;
TitleColor = 'k';

clr.paper = 1;
clr.title = 'k';
clr.bck = 'w';
clr.ax = 'k';

% savingDir = 'D:\ArthurC\Fight_Club\PARLA_figures';
myfilepath_res = [workingdir filesep 'img_paper4' filesep];
if exist(myfilepath_res)~=7
    mkdir(myfilepath_res)
end

% savingDir = [workingdir filesep 'img_paper_def'];


%% FIG 5: MATOUT target %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig39 = figure(39);clf
fig39.Position = [20 150 size(MatOutTarget,2)/2 size(MatOutTarget,1)/2];
% fig40.Position(4) = fig40.Position(3)*fig40.Position(4)/1883;
fig39.InvertHardcopy = 'off';
fig39.Color = clr.bck;

tight_subplot(1,1,0,0,0.0);

im = imagesc(imgaussfilt((MatOutTarget).^(1/2),1));
% tname = text(12, 12,'Target MatOut','color','w','fontsize',labelSize-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
axis image, colormap hot, axis off

posScale(1) = 60; posScale(2) = size(MatOutTarget,1)-70;
hold on
pscale = plot(posScale(1)+[1 100],posScale(2).*[1 1],'w-','linewidth',3);
tscale = text(mean(pscale.XData), pscale.YData(1),'100µm','color','w','fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center');
im.CData(posScale(2)*[1 1],posScale(1)+[1:100])=mean(caxis);
    

figname = 'fig5_MatOutTarget';
saveas(fig39,[myfilepath_res filesep figname '.png'])
savefig(fig39,[myfilepath_res filesep figname])
WriteTif(im.CData,hot,[myfilepath_res filesep figname '.tif'])

%% FIG 5: MATOUTS zoom horseshoe %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([myfilepathdata '_MatOut_multi_' num2str(abs(30)) 'dB'],'MatOut','ULM','P','listAlgo','CountMatOut');
%%
fig40 = figure(40);clf
fig40.Position = [169 511 688 426];
fig40.InvertHardcopy = 'off';fig40.Color = clr.bck;
a = tight_subplot(2,3,[0.02 0.005],0,0.0);

list2Disp = [1 2 3 6 7 8];
xcrop = [230:360];
ycrop = [210:330];
figname = 'fig5_MatOut_horseshoe_30dB';
for ii=1:numel(list2Disp)
    ialgo = list2Disp(ii);
    axes(a(ii))
    if ii ==1
        im = imagesc(imgaussfilt((MatOutTarget(ycrop,xcrop)).^(1/2),.001));
        tname(ialgo) = text(3, 3,'Target MatOut','color','w','fontsize',labelSize-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    else
        im = imagesc(imgaussfilt((MatOut{ialgo-1}(ycrop,xcrop)).^(1/2),.001));
        tname(ialgo) = text(3, 3,listAlgoName{ialgo-1},'color','w','fontsize',labelSize-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    end
    axis image, colormap hot
    if ii==1
        posScale(1) = 8; posScale(2) = size(im.CData,1)-5;
        hold on
        im.CData(posScale(2)*[1 1],posScale(1)+[1:50])=mean(caxis);
        tscale = text(posScale(1)+25, posScale(2),'50µm','color','w','fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center');
    end
    
    a0 = gca;a0.XColor = clr.ax;a0.YColor = clr.ax;a0.XTick = []; a0.YTick = [];
    if ii==1
        WriteTif(im.CData,hot,[myfilepath_res filesep figname 'Target' '.tif'])
    else
        WriteTif(im.CData,hot,[myfilepath_res filesep figname ShortName{ialgo-1} '.tif'])
    end
end

% Save fig
figname = 'fig5_MatOut_horseshoe_30dB';
saveas(fig40,[myfilepath_res filesep figname '.png'])
savefig(fig40,[myfilepath_res filesep figname])

%% FIG 5: MATOUTS zoom 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([myfilepathdata '_MatOut_multi_' num2str(abs(30)) 'dB'],'MatOut','ULM','P','listAlgo','CountMatOut');
%%
fig42 = figure(41);clf
fig42.Position = [169 257 1133 682];
fig42.InvertHardcopy = 'off';fig42.Color = clr.bck;
a = tight_subplot(2,3,[0.01 0.01],0,0.0);

% list2Disp = [1 2 6 7];
xcrop = [250:650];
ycrop = [350:590];

figname = 'fig5_MatOut_zoom_30dB';
for ii=1:numel(list2Disp)
    ialgo = list2Disp(ii);
    axes(a(ii))
    if ii ==1
        im = imagesc(imgaussfilt((MatOutTarget(ycrop,xcrop)).^(1/2),.5));
        tname(ialgo) = text(5, 5,'Target MatOut','color','w','fontsize',labelSize+2-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    else
        im = imagesc(imgaussfilt((MatOut{ialgo-1}(ycrop,xcrop)).^(1/2),.001));
        tname(ialgo) = text(5, 5,listAlgoName{ialgo-1},'color','w','fontsize',labelSize+2-2,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    end
    axis image, colormap hot

    if ii==1
        posScale(1) = 10; posScale(2) = size(im.CData,1)-5;
        im.CData(posScale(2)*[1 1],posScale(1)+[1:50])=mean(caxis);

        hold on
        pscale = plot(posScale(1)+[1 50],posScale(2).*[1 1],'w-','linewidth',3);
        tscale = text(mean(pscale.XData), pscale.YData(1),'50µm','color','w','fontsize',labelSize+2-2,'VerticalAlignment','bot','HorizontalAlignment','center');
    end
    
    a0 = gca;a0.XColor = clr.ax;a0.YColor = clr.ax;a0.XTick = []; a0.YTick = [];
    if ii==1
        WriteTif(im.CData,hot,[myfilepath_res filesep figname 'Target' '.tif'])
    else
        WriteTif(im.CData,hot,[myfilepath_res filesep figname ShortName{ialgo-1} '.tif'])
    end
end

% Save fig
figname = 'fig5_MatOut_zoom_30dB';
% saveas(fig42,[myfilepath_res filesep figname '.png'])
% savefig(fig42,[myfilepath_res filesep figname])

%% FIG 3: Error lat ax flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([myfilepathdata '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')

% remove biais
if 0
    for ialgo=1:numel(listAlgo)
        ErrList{ialgo}(:,2)=ErrList{ialgo}(:,2)-mean(ErrList{ialgo}(:,2));
        ErrList{ialgo}(:,3)=ErrList{ialgo}(:,3)-mean(ErrList{ialgo}(:,3));
        ErrList{ialgo}(:,1)=vecnorm([ErrList{ialgo}(:,3),ErrList{ialgo}(:,2)],2,2);
    end
end
%%
fig43 = figure(20);clf
fig43.Position = [223 414 1250 340];
fig43.InvertHardcopy = 'off';fig43.Color = clr.bck;
postxt_x = -.4;

aa0 = tight_subplot(2,7,[.09 .015],[.07 .01],[.01 .01]);
xedge = linspace(-.3,.3,100);

labelSize = 16;
for ii = 1:numel(listAlgo)
    ialgo=ii;
    axes(aa0(ii))
    
    axhist_x = histogram((ErrList{ialgo}(:,3)),xedge,'edgecolor','none','Normalization','Probability','FaceColor',[1 1 1]*0.2);
    axhist_x.FaceAlpha = 1;
    axhist_x.FaceColor = ListColor(ii,:);
    ErrMean_x(ialgo) = mean(ErrList{ialgo}(:,3));
    ErrStd_x(ialgo) = std(ErrList{ialgo}(:,3));

    hold on
    xlim([-1 1]*.3),ylim([0 .06])
    
    text(.02,aa0(ii).YLim(2)-.009,...
        ['\sigma=' num2str(round(ErrStd_x(ialgo),2),'%.2f') ],...
        'fontsize',labelSize,'HorizontalAlignment','left','VerticalAlignment','top','color',clr.title)
    
    grid on
    
    ax = gca;
    ax.Color = 'none';
    ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';
    ax.YTickLabel = [];
    ax.XTickLabel = {'-\lambda/4','0','\lambda/4'};ax.XTick = [-.25,0,.25];
    ax.FontSize = labelSize-4;
    ax.TickDir = 'in';
    ax.TickLength = [.04 .2];
    ax.Box='off';
    ax.GridAlpha = .5; ax.GridLineStyle = '--';
%     title(listAlgoName{ialgo},'color',clr.title,'fontsize',labelSize,'interpreter','none')
    if ii == 1
        text(postxt_x,mean(ylim),'Lateral error','color',clr.title,'FontSize',labelSize,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90)
    end
end

%
for ii= 1:numel(listAlgo)
    ialgo=ii;
    axes(aa0(ii+7))

    axhist_z = histogram(ErrList{ialgo}(:,2),xedge,'edgecolor','none','Normalization','Probability','FaceColor',[1 1 1]*0.2);
    axhist_z.FaceAlpha = axhist_x.FaceAlpha;
    ErrMean_z(ialgo) = mean(ErrList{ialgo}(:,2));
    ErrStd_z(ialgo) = std(ErrList{ialgo}(:,2));
    axhist_z.FaceColor = ListColor(ii,:);

    hold on
    xlim([-1 1]*.3)
    ylim([0 .03])
    
    text(.02,aa0(ii+7).YLim(2)-.001,...
        ['\sigma=' num2str(round(ErrStd_z(ialgo),2),'%.2f') ],...
        'fontsize',labelSize,'HorizontalAlignment','left','VerticalAlignment','top','color',clr.title)    
    
    grid on
    set(gca, 'yticklabel', []);set(gca, 'xticklabel', {'-\lambda','0','+\lambda'});
    
    az = gca;
    az.Color = ax.Color;
    az.XAxisLocation = 'origin';az.YAxisLocation = 'origin';
    az.YTickLabel = [];
    az.XTickLabel = ax.XTickLabel;az.XTick = ax.XTick;
    az.FontSize = ax.FontSize;
    az.TickDir = ax.TickDir;
    az.TickLength = ax.TickLength;
    az.Box= ax.Box;
    az.GridAlpha = ax.GridAlpha; az.GridLineStyle = ax.GridLineStyle;
    if ii ==1
        text(postxt_x,mean(ylim),'Axial error','color',clr.title,'FontSize',labelSize,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90)
    end
end

% Save fig
figname = 'fig3_distrib_30dB';
% savefig(fig43,[myfilepath_res filesep figname])
print(fig43,[myfilepath_res filesep figname],'-dpng','-r600')

%%
ErrMean_x_interp = mean([ErrList{3}(:,3);ErrList{4}(:,3);ErrList{5}(:,3)]);
ErrStd_x_interp = std([ErrList{3}(:,3);ErrList{4}(:,3);ErrList{5}(:,3)]);
ErrMean_z_interp = mean([ErrList{3}(:,2);ErrList{4}(:,2);ErrList{5}(:,2)]);
ErrStd_z_interp = std([ErrList{3}(:,2);ErrList{4}(:,2);ErrList{5}(:,2)]);


%% FIG 3: TRUE POSITIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([myfilepathdata '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')

for ialgo=1:numel(listAlgo)
TP(ialgo) = Stat_class(3,ialgo);
FN(ialgo) = Stat_class(4,ialgo);
FP(ialgo) = Stat_class(5,ialgo);

end

J_p = TP./(FP+TP)*100; % precision
J_r = TP./(FN+TP)*100; % sensitivity
J_ac = TP./(FP+FN+TP)*100; % Jaccard

%%
fig50=figure(50);clf;hold on
fig50.Position = [880 600 708 250];
fig50.InvertHardcopy = 'off';fig50.Color = clr.bck;
b1(1) = barh(1:7, [TP']*1e-3, 1, 'stack','DisplayName','True Positive');
b1(2:3) = barh(1:7, [ -FN', -FP']*1e-3, 1, 'stack');

b1(1).FaceColor = 'flat';
b1(1).CData = ListColor;
b1(1).FaceAlpha =1;
b1(2).FaceColor = 'flat';
b1(2).CData = ListColor;
b1(2).FaceAlpha =.5;
b1(3).FaceColor = 'flat';
b1(3).CData = ListColor;
b1(3).FaceAlpha =.2;

b1(2).DisplayName = 'False Negative';b1(3).DisplayName = 'False Positive';
legend('box','off','Location','northwest')

for ii=1:numel(listAlgo)
   text(TP(ii)*1e-3+2e1,ii,listAlgoName(ii))  
   text(TP(ii)*1e-3/2,ii, [num2str(round(TP(ii)/1000)) 'k'],'HorizontalAlignment','center')
   text(-FN(ii)*1e-3/2,ii, [num2str(round(FN(ii)/1000)) 'k'],'HorizontalAlignment','center')
   text(-FP(ii)*1e-3/2-FN(ii)*1e-3,ii, [num2str(round(FP(ii)/1000)) 'k'],'HorizontalAlignment','center')
end
xlim([-7*1e2 5*1e2])
ylim([.5 7.5])

a =gca;
a.YTick =[];
% set(gca,'ytick',[])
a.YAxisLocation = 'origin';

a.XAxis.Exponent = 0;
xtickformat('%.0fk')
a.XAxis.TickDirection = 'out';
a.Position = [.05 .12 .9 .87];
% % % ytickformat('%g%%')

text(a.XLim(1),.5,'False','HorizontalAlignment','right','VerticalAlignment','middle')
text(a.XLim(2),.5,'True','HorizontalAlignment','left','VerticalAlignment','middle')

NbDetectableBulles = mean(ValTP(1,:)+ValTP(2,:));
plot([1 1]*NbDetectableBulles/1000,[0 7.5],'--','LineWidth',1,'color',[0 0 0 .2],'HandleVisibility','off')
text(NbDetectableBulles/1000,.2,'Max','HorizontalAlignment','center')
grid on

%% Save fig
figname = 'fig3_TP_FP_FN_30dB';
print(fig50,[myfilepath_res filesep figname],'-dpng','-r600')
% savefig(fig50,[myfilepath_res filesep figname])

%% FIG 3: JACCARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig52=figure(52);clf
fig52.Position = [223 590 875 260];
fig52.InvertHardcopy = 'off';fig52.Color = clr.bck;
barcolor = [1 1 1]*.5;

stsuba= tight_subplot(1,3,[0 .02],[.15 .15],[.01 0.05]);

axes(stsuba(1));

b_p = barh(J_p,'facecolor',barcolor,'BarWidth',1);
text(b_p.XData,b_p.YData+20,string(round(b_p.YData,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top','color',clr.title,'FontSize',labelSize)
b_p.FaceColor = 'flat';
b_p.CData = ListColor;
text(J_p-1,1:length(J_p),num2str(round(J_p,1)'),'VerticalAlignment','middle','HorizontalAlignment','right'); 

ylim([0.5 7.5])
xlim([0 90])
ap = gca;
ap.XAxisLocation = 'origin'; ap.YAxisLocation = 'origin';
ap.YTickLabel = {};
ap.FontSize = labelSize-2;
ap.TickDir = 'in';
ap.Box = 'off';
ap.XAxis.FontSize = labelSize -4;grid on

text(5,7.5,'Precision','FontSize',labelSize,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold')

axes(stsuba(2));suba = gca;

b_r = barh(J_r,'facecolor',barcolor,'BarWidth',1);
text(b_r.XData,b_r.YData+20,string(round(b_r.YData,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top','color',clr.title,'FontSize',labelSize)
b_r.FaceColor = 'flat';
b_r.CData = ListColor;
text(J_r-1,1:length(J_r),num2str(round(J_r,1)'),'VerticalAlignment','middle','HorizontalAlignment','right'); 


ylim([0.5 7.5])
xlim([0 80])

ar = gca;
ar.XAxisLocation = 'origin';ar.YAxisLocation = 'origin';

ar.YTickLabel = {};
ar.FontSize = ap.FontSize;
ar.TickDir = ap.TickDir;
ar.XAxis.FontSize = ap.XAxis.FontSize;
ar.Box= 'off';grid on

text(5,7.5,'Sensitivity','FontSize',labelSize,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold')

axes(stsuba(3));suba = gca;

b_ac = barh(J_ac,'facecolor',barcolor,'BarWidth',1);
text(b_ac.XData,b_ac.YData+20,string(round(b_ac.YData,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top','color',clr.title,'FontSize',labelSize)
b_ac.FaceColor = 'flat';
b_ac.CData = ListColor;
text(J_ac-1,1:length(J_ac),num2str(round(J_ac,1)'),'VerticalAlignment','middle','HorizontalAlignment','right'); 

ylim([0.5 7.5])
xlim([0 100])

aa = gca;
aa.XAxisLocation = 'origin';aa.YAxisLocation = 'origin';
aa.FontSize = ap.FontSize;
aa.TickDir = ap.TickDir;
aa.XAxis.FontSize = ap.XAxis.FontSize;
aa.YTickLabel = {};
aa.Box = 'off';grid on

text(b_ac.YData+5,b_ac.XData,listAlgoName(1:7),...
    'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize-2)

text(5,7.5,'Jaccard Index','FontSize',labelSize,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold')

%% Save fig
figname = 'fig3_Jacc_30dB';
print(fig52,[myfilepath_res filesep figname],'-dpng','-r600')
% savefig(fig52,[myfilepath_res filesep figname])

%% %%%%%%% FIG 2 : BRANCHES MATOUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MM = MatOutTarget(780:800,120:800);

DistNul = [];Gap = [];
ListGap = [];
seuil_p=.03;
% seuil=5;
seuil = seuil_p.*max(MM(:,1));
Nresize = 1;
figure(1)
for ialgo = 1:8
    ialgo
    if ialgo <8
        MM = MatOut{ialgo}(780:800,220:800);
    else
        MM = MatOutTarget(780:800,220:800);
    end    
    lz = 1:2:size(MM,2);
    for ind=1:numel(lz)
        Lseq = MM(:,lz(ind));
        Lseq(Lseq<0)=0;
        Lseq = Lseq(find(Lseq>seuil,1,'first'):find(Lseq>seuil,1,'last'));
        DistNul(ind,ialgo) = sum(Lseq<seuil);
        
        if ialgo==8
            im = imregionalmax(Lseq);
            posMax = find(im); Kepmax = [0;diff(posMax)]~=1; posMax=posMax(Kepmax);
            val = Lseq(im);
            [b1,b2] = sort(val(Kepmax),'descend');
            try
                Gap(ind) = abs(posMax(b2(1))-posMax(b2(2)));
            catch
                Gap(ind)=0;
            end
        end
    end
    
    [mingap,~] = find(DistNul(:,ialgo)==0,1,'last');
    ListGap(ialgo) = mingap+1;
end
SR2conv = 1./ULM.res; % convert in \lambda

%% explanation figure
fig27 = figure(27);clf
ialgo = 7; zpos = 350;
MM = MatOut{ialgo}(780:800,220:800);
subplot 211;hold on
imagesc(MM.^(1/2));axis off; colormap hot;colorbar
Lseq = MM(:,zpos);

thres = max(MM(:,1))*seuil_p;
% seuil = thres+1+1;

l_x = 1:numel(Lseq);
Lseqp = imresize(Lseq,[numel(Lseq)*Nresize,1]);
plot(zpos*[1 1],[1 size(MM,1)],'w:','LineWidth',2)
title(['MatOut: ' listAlgoName{ialgo}])

subplot 212;hold on
plot(l_x*SR2conv,Lseq,'-','MarkerSize',7,'LineWidth',1.1,'Marker','.')
plot([l_x(1) l_x(end)]*SR2conv,[1 1]*seuil,'r--','LineWidth',1)

Pic(1) = find(Lseq>seuil,1,'first');
Pic(2) = find(Lseq>seuil,1,'last');
RangeGap = find(Lseq(Pic(1):Pic(2))<seuil);
pp = plot(l_x(RangeGap([1 end]) + Pic(1)-1)*SR2conv,[1 1]*seuil,'r','LineWidth',3);

text(mean(pp.XData)*SR2conv,mean(pp.YData)+1,'Gap','VerticalAlignment','bot','HorizontalAlignment','center','FontSize',labelSize)
grid on
legend({['Interpolated profil (x' num2str(Nresize) ')'],'MatOut profil',['Threshold (' num2str(seuil) ')']})
xlim((l_x(Pic) + [-1 2]*2)*SR2conv)
xlabel('[\lambda]');ylabel('Intensity')
title('Profil')
suptitle('Estimation of the Gap.')

% saveas(fig27,[savingDir filesep filename '_ExplainGap.png'])
%%
DistNul_wv = DistNul*SR2conv;
fig26 = figure(26);
fig26.Position = [500 200 1052 493];
clf; clear plt_dist
ax1 = axes('Position',[.05 .35 .92 .6]);

Gap = DistNul_wv(:,end); % take the measured gap on the Traget MatOut in [nbSRpixel]
pp = polyfit(lz(:),Gap(:),1);

Gap_pente = pp(2);
Gap_origine = pp(1);
Gap_linear = pp(2) + lz*pp(1);

hold on
% plDist = plot([0 Gap_linear(end)],[0 Gap_linear(end)],'k:','LineWidth',2,'DisplayName','Distance between\newlinecanal centers');

ListGap_pix = Gap_linear(ListGap);

jumpMarker = 3;
for ialgo=1:7
    plt_dist(ialgo) = plot(Gap_linear(:),DistNul_wv(:,ialgo),'Marker',ListMarker{ialgo},'MarkerSize',7,'DisplayName',listAlgoName{ialgo},'color',ListColor(ialgo,:),'linewidth',1.5);
    plt_dist(ialgo).MarkerFaceColor = plt_dist(ialgo).Color;
    plt_dist(ialgo).MarkerIndices = plt_dist(ialgo).MarkerIndices((ialgo*3):20:end);
    
    plot(ListGap_pix(ialgo),0,'+','color',ListColor(ialgo,:),'MarkerSize',10)
    text(ListGap_pix(ialgo),0,num2str(round(ListGap_pix(ialgo),2)),'VerticalAlignment','top','HorizontalAlignment','center')
end

lgd = legend([plt_dist],'Location','northwest','Interpreter','tex');
% lgd.Box = 'off';
xlim([3 9]*SR2conv)
ylim([0 12]*SR2conv)
ylabel('Measured canal to canal distance [\lambda]')
xlabel('Simulated canal to canal distance')
grid on 
title(['Size of the measured gap on Matout (threshold ' num2str(seuil) ', ' num2str(seuil_p*100) '%)'])

for ialgo=0:7
    ax(ialgo+1) = axes('Position',[.05 + (ialgo)/7*.8,.01,.8/8,.2]);
    if ialgo >0
        MM = MatOut{ialgo}(780:800,220:800);
    else
        MM = MatOutTarget(780:800,220:800);
    end
    imagesc(MM.^(1/2)),axis off, colormap hot
    if ialgo>0,title(listAlgoName{ialgo},'Interpreter','none');else title('Target','Interpreter','none');end
end

%% Save fig
figname = 'fig2_gap_30dB';
print(fig26,[myfilepath_res filesep figname],'-dpng','-r600')
savefig(fig26,[myfilepath_res filesep figname])

%%

%% explanation figure
fig28 = figure(28);clf
fig28.Position = [290 439 680 500];
axes('Position',[.0 .5 1 .5])
ialgo = 7;
MM = MatOut{ialgo}(780:800,220:800);
hold on
MM = MM.^(1/2);
nslice = size(MM,1);
lx = ones(1,nslice);ly = (1:nslice)*SR2conv;
for zz=30:40:size(MM,2)   
    plot3(lx*zz*SR2conv,ly, MM(:,zz),'k-','LineWidth',1.5)
%     text(zz,1,1,num2str(zz),'VerticalAlignment','bot','FontSize',8)
end

zz0 = 350;
plot3(lx*zz0*SR2conv,ly, MM(:,zz0),'k-','LineWidth',3)
text(zz0,1,1,num2str(zz0),'VerticalAlignment','bot','FontSize',10)


% ss = surf(MM,'AlphaData',MM,'EdgeColor','none','FaceAlpha','flat');
% ss.AlphaData =MM/max(MM(:));
% ss.FaceAlpha = 'interp'
grid on
myclr = flip(hot,1);
imagesc([1:size(MM,2)]*SR2conv,[1:size(MM,1)]*SR2conv,imgaussfilt(MM,.5));
colormap(myclr),colorbar
view([-144 22])
pbaspect([3 1 1/3])
title(listAlgoName{ialgo})
ylabel('[\lambda]'),xlabel('[\lambda]')
axes('Position',[.1 .1 .8 .3])

hold on
Lseq = MM(:,zpos);

thres = max(MM(:,1))*seuil_p;
% seuil = thres+1+1;

l_x = (1:numel(Lseq))*SR2conv;
plot(zpos*[1 1],[1 size(MM,1)],'w:','LineWidth',2)
title(['MatOut: ' listAlgoName{ialgo}])

plot(l_x,Lseq,'LineWidth',1.1)
plot(l_x,Lseq,'k.','MarkerSize',7)
plot([l_x(1) l_x(end)],[1 1]*seuil,'r--','LineWidth',1)

Pic(1) = find(Lseq>seuil,1,'first');
Pic(2) = find(Lseq>seuil,1,'last');
RangeGap = find(Lseq(Pic(1):Pic(2))<seuil);
pp = plot(l_x(RangeGap([1 end]) + Pic(1)-1),[1 1]*seuil,'r','LineWidth',3);

text(mean(pp.XData),mean(pp.YData)+1,'Gap','VerticalAlignment','bot','HorizontalAlignment','center','FontSize',labelSize)
grid on
legend({['Interpolated profil (x' num2str(Nresize) ')'],'MatOut profil',['Threshold (' num2str(seuil) ')']})
xlim(l_x(Pic) + [-1 2]*2*SR2conv)
xlabel('[\lambda]');ylabel('Intensity')
title('Profil')
suptitle('Estimation of the Gap.')

%%
saveas(fig28,[myfilepath_res filesep filename '_gapSlice.png'])
savefig(fig28,[myfilepath_res filesep filename '_gapSlice'])


%% Fig 2 : RMSE(dB) silico
load([myfilepathdata '_Stats_multi' num2str(abs(30)) 'dB' '.mat'],'-mat')

% remove specific error
val_Err = {};RMSE_m=[];
for ii=1:numel(ClutterList) 
    ii
    temp = load([myfilepathdata '_Stats_multi' num2str(abs(ClutterList(ii))) 'dB' '.mat'],'-mat');
    for ialgo=1:numel(listAlgo)
        temp.ErrList{ialgo}(:,2)=temp.ErrList{ialgo}(:,2)-mean(temp.ErrList{ialgo}(:,2));
        temp.ErrList{ialgo}(:,3)=temp.ErrList{ialgo}(:,3)-mean(temp.ErrList{ialgo}(:,3));
        temp.ErrList{ialgo}(:,1)=vecnorm([temp.ErrList{ialgo}(:,3),temp.ErrList{ialgo}(:,2)],2,2);
        val_Err{ialgo,ii} = vecnorm([temp.ErrList{ialgo}(:,3),temp.ErrList{ialgo}(:,2)],2,2);
        
        if ClutterList(ii)==-30
            RMSE_m(ialgo) = mean(val_Err{ialgo,ii});
        end 
    end
end

%%
fig60 = figure(60);clf
fig60.Position = [190 470 1150 240];
axes('Position',[.04 .05 .96 .95])
vv = myBoxPlotFC_multi(listAlgoName,val_Err,15,ListColor',10,0);
ylim([0 .26])

for ii=1:numel(ClutterList)
    tt=text(.55+(ii-1)/(numel(ClutterList)-1)*.9,.01,[num2str(ClutterList(ii))],'VerticalAlignment','bot','HorizontalAlignment','center','FontSize',9);
end
tt=text(.55+(ii)/(numel(ClutterList)-1)*.9,.01,'dB','VerticalAlignment','bot','HorizontalAlignment','center','FontSize',9);

% ylabel('RMSE')
grid on
a = gca;a.FontSize = 11;a.YTick = [0:.05:.25];
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;

%% Save fig
figname = 'fig2_RMSEsilico';
% savefig(fig60,[myfilepath_res filesep figname])
print(fig60,[myfilepath_res filesep figname],'-dpng','-r600')

%% Save data for chart

Radar.RMSE = RMSE_m;    
Radar.Jaccard = J_ac;
Radar.precision = J_p;
Radar.Gap       = ListGap_pix;
Radar.listName = listAlgoName;
save([workingdir filesep filename '_radar_simu'],'Radar')


