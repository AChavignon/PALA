%% ULM PARLA: Post Processing - displaying and anylisng results for IN VIVO BRAIN
% For each algorithm, bubbles are detected, localised and tracks.
% This script compare results of differents algorithms building matout, and comparing various behaviours. 
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
PARLA_data_folder = 'F:\ArthurC\Data\Fight_Club\Data'; % path of data
addpath(genpath(PARLA_addons_folder))
%% Selected data file and saving folders
workingdir = [PARLA_data_folder '\PARLA_data_InVivoBrain'];
filename = 'PARLA_InVivoRatBrain_';
cd(workingdir)

SavingFolder = 'fightclub_round2';
myfilepath = [workingdir filesep filename];
savingpath = [workingdir filesep SavingFolder filesep filename];

irecon = '1';
load([savingpath 'MatOut_multi' irecon '.mat'])
load([savingpath 'MatOut_multi_nointerp' irecon '.mat'],'ProcessingTime','ULM','listAlgo','MatOutNoInterp')
[listAlgoName,ListColor,ListMarker,ListShortName] = GetFormatFightClub;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelSize = 15;
TitleColor = 'k';

clr.paper = 1;
clr.title = 'k';
clr.bck = 'w';
clr.ax = 'k';
fol = '';

if exist([workingdir filesep SavingFolder filesep 'img_paper4' fol])~=7
    mkdir([workingdir filesep SavingFolder filesep 'img_paper4' fol])
end
widthBar = .4;

savingDir = [workingdir filesep SavingFolder filesep 'img_paper4' fol];

%% FIGURE MATOUT + BULLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig39 = figure(39);clf
fig39.Position = [750 230 1080 700];
fig39.InvertHardcopy = 'off';
fig39.Color = clr.bck;

hold on
% load IQ
if ~exist('bulles')
    load([workingdir filesep 'IQ' filesep filename '_001.mat'],'IQ')
    bulles = SVDfilter(IQ,ULM.seuil);
    % bulles = filter(but_b,but_a,bulles,[],3);
    bulles(~isfinite(bulles))=0;
end
dB = 20*log10(abs(bulles(1:end,1:end-10,:))); dB = dB-max(dB(:));
l_z = 1:size(dB,1);l_x = 1:size(dB,2);
imagesc(l_x,l_z,dB(:,:,20),[-40 0]),colormap gray
axis image;axis off

Mat_tmp =imgaussfilt(MatOut{end}(17:765,23:1166).^(1/3),.01);

h=colorbar;
colorTitleHandle = get(h,'Title');set(colorTitleHandle ,'String','dB');

Mat_tmp = Mat_tmp./max(Mat_tmp(:));
Mat_rgb = ind2rgb(round(Mat_tmp*63),hot);

im = image(l_x,l_z,Mat_rgb);

AlphaMap = zeros(size(Mat_rgb(:,:,1)));
ex = (1:size(AlphaMap,2))-size(AlphaMap,2)/2-100;
ez = (1:size(AlphaMap,1))-size(AlphaMap,1)/2;
[ex,ez]=meshgrid(ex,ez);

AlphaMap = -ex*3 - ez + size(AlphaMap,2)/3;
AlphaMap = AlphaMap./600;
AlphaMap(AlphaMap>1)=1;AlphaMap(AlphaMap<0)=0;

im.AlphaData = AlphaMap*.9;

posScale(1) = 2; posScale(2) = 4;
hold on
pscale = plot(posScale(1)+[1 10],posScale(2).*[1 1],'w-','linewidth',3);
tscale = text(mean(pscale.XData), pscale.YData(1),'1mm','color',clr.bck,'fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center');

ca = gca;
ca.YDir = 'reverse';
ca.Position=[.005 .01 .94 .96];

saveas(fig39,[savingDir filesep filename '_MatBulles' irecon '.png'])
savefig(fig39,[savingDir filesep filename '_MatBulles' irecon ])

%% %% FIG 6 : MATOUT RAIDAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig38 = figure(40);clf
fig38.Position = [5 160 1900 830];
fig38.InvertHardcopy = 'off';
fig38.Color = clr.bck;
countpos = [-20 -20];

ialgo = 7;

ii=imagesc(imgaussfilt(MatOut{ialgo}(17:765,23:1166).^(1/3),.01));

fig38.Position(4:-1:3)=size(ii.CData);
a = gca;
a.Position = [0 0 1 1];
axis image, colormap hot
text(12,10,listAlgoName{ialgo},'color','w','fontsize',labelSize,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
axis off
    
posScale(1) = 60; posScale(2) = size(MatOut{ialgo},1)-70;
hold on
pscale = plot(posScale(1)+[1 100],posScale(2).*[1 1],'w-','linewidth',3);
tscale = text(mean(pscale.XData), pscale.YData(1),'1mm','color',clr.bck,'fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center');
countpos = size(ii.CData)+countpos;
ii.CData(posScale(2)*[1 1],posScale(1)+[1:100])=mean(caxis);

caxis([0 10])
tcount = text(countpos(2), countpos(1),[num2str(round(NbrOfLoc(ialgo)/1e6,2)),'M'],'color',clr.bck,'fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','right');

figname = 'fig6_matout_radial_full';
% saveas(fig38,[savingDir filesep figname '.png'])

aa = gca;
WriteTif(ii.CData,aa.Colormap,[savingDir filesep figname '.tif'],'caxis',caxis)

%% %% FIG 6 : MATOUTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig40 = figure(40);clf
fig40.Position = [5 160 1290 860];
fig40.InvertHardcopy = 'off';
fig40.Color = clr.bck;

% a = tight_subplot(2,4,[0.00 0.02],0.02,0.01);
a = tight_subplot(2,2,[0.01 0.005],0.00,0.00);
countpos = [-20 -20];
listAlgoMatOut = [2 3 5 7];
extradata= 20*0;
for ii=1:numel(listAlgoMatOut)
    ialgo=listAlgoMatOut(ii);
    
    axes(a(ii))
    im=imagesc((MatOut{ialgo}((-extradata+476):(627+extradata),(-extradata+651):(882+extradata))).^(1/3));
    axis image, colormap hot
    
    axis off
    
    if ii==1
        posScale(1) = 10+extradata;
        posScale(2) = round(a(ii).YLim(2)-10);
        hold on
        pscale = plot(posScale(1)+[1 50],posScale(2).*[1 1],'w-','linewidth',3);
        tscale = text(mean(pscale.XData), pscale.YData(1),'500µm','color',clr.bck,'fontsize',labelSize+4,'VerticalAlignment','bot','HorizontalAlignment','center');
        countpos = [a(ii).YLim(2) a(ii).XLim(1)]-[5,10];
        countpos = [a(ii).YLim(2) a(ii).XLim(2)]+[-5,-10];
        im.CData(posScale(2)*[1 1],posScale(1)+[1:50])=mean(caxis);
    end
    caxis([0 10])
    text(extradata+2,extradata+2,listAlgoName{ialgo},'color','w','fontsize',labelSize+4,'VerticalAlignment','top','HorizontalAlignment','left','backgroundcolor','k');
    tcount = text(countpos(2), countpos(1),[num2str(round(NbrOfLoc(ialgo)/1e6,2)),'M'],'color',clr.bck,'fontsize',labelSize+4-2,'VerticalAlignment','bot','HorizontalAlignment','right');
    
    WriteTif(im.CData,hot,[savingDir filesep figname num2str(ii) '.tif'],'caxis',caxis)
end
axes(a(end)),axis off
linkaxes(a)

%% Save fig
figname = 'fig6_matouts_tight';
saveas(fig40,[savingDir filesep figname '.png'])
savefig(fig40,[savingDir filesep figname])

%% %% FIG 6 : ALIASING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Aliasing Index with MatOut WITHOUT interpolation

T_alias = 10;
F_alias = 1/T_alias;
harm = [1 2]';
Freq_i_x = 119*harm + [-1:1];
Freq_base_x = 119*harm + [-5:-3,3:5];
Freq_i_z = 79*harm + [-1:1];
Freq_base_z = 79*harm + [-5:-3,3:5];
clear TF2_z_list TF2_x_list
for ialgo=1:numel(listAlgo)

    TF2 = (fft2(MatOutNoInterp{ialgo}));
%     TF2 = (fft2(MatOut{ialgo}));

    TF2_x = smooth(sum(abs(TF2).^2,1)/sum(sum(abs(TF2).^2,1)),1);
    TF2_z = smooth(sum(abs(TF2).^2,2)/sum(sum(abs(TF2).^2,1)),1);
    
%     TF2_z = TF2_z(1:(size(TF2_z,1)-1)/2) + flip( TF2_z((size(TF2_z,1)-1)/2+2:end) ); 
%     TF2_x = TF2_x(1:(size(TF2_x,1)-1)/2) + flip( TF2_x((size(TF2_x,1)-1)/2+2:end) ); 
    
    
    peak_z = max(TF2_z(Freq_i_z),[],2);
    baseline_z = mean(TF2_z(Freq_base_z),2);
    P_alias_z(ialgo) = 20*log10(sum(peak_z./baseline_z));
    
    peak_x = max(TF2_x(Freq_i_x),[],2);
    baseline_x = mean(TF2_x(Freq_base_x),2);
    P_alias_x(ialgo) = 20*log10(sum(peak_x./baseline_x));
    
    P_alias(ialgo) = 20*log10(1/2*sum(peak_x./baseline_x  + peak_z./baseline_z));
    
    TF2_z_list(ialgo,:) = TF2_z;
    TF2_x_list(ialgo,:) = TF2_x;
 
end

%% Build figure plot

figure(71);clf
subplot 211
hold on
for ialgo=1:numel(listAlgo)
    plot(20*log10(TF2_z_list(ialgo,:)))
end
legend(listAlgoName)
subplot 212
hold on
for ialgo=1:numel(listAlgo)
    plot(20*log10(TF2_x_list(ialgo,:)))
end
legend(listAlgoName)

%% Build figure
clear ae
fig41 = figure(41);clf
fig41.Position = [530 350 670 260];
fig41.InvertHardcopy = 'off';
fig41.Color = clr.bck;

P_alias_bis = P_alias;
% P_alias_bis(2) = 20;


b_al = barh(P_alias_bis,'facecolor',[1 1 1]*0.5,'BarWidth',1);
% text(b_al.XData,b_al.YData,string(round(b_al.YData,1)),...
%     'HorizontalAlignment','center','VerticalAlignment','bot','color','w','FontSize',labelSize)

b_al.FaceColor = 'flat';
b_al.CData = ListColor;

text(b_al.YData-.1,b_al.XData,string(round(P_alias,1)),...
    'HorizontalAlignment','right','VerticalAlignment','middle','color',clr.title,'FontSize',labelSize-2)

text(b_al.YData+.5,b_al.XData,listAlgoName,...
    'HorizontalAlignment','left','VerticalAlignment','middle','color',clr.title,'FontSize',labelSize)


ae = gcf;ae = ae.Children;
% ae.YLim = [0 18];
% grid on
ylim([0 8] + [1 -1]*.5)
xlim([0 23])
hold on
% plot([-1 1]*.4+ 2,[0 1]+P_alias_bis(2)-2,clr.bck,'LineWidth',8)

ae.Position = [.01 .11 .9 .78];
ae.Color = clr.bck;

ae.XAxisLocation = 'origin'; ae.YAxisLocation = 'origin';
ae.XAxis.Color = clr.ax;
ae.YAxis.Color = ae.XAxis.Color;
% ae.XTickLabel = listAlgoName2;
ae.FontSize = labelSize-2;
ae.XAxis.FontSize = labelSize;
ae.TickDir = 'both';
ae.TickLength = [.005 .0];
ae.Box= 'off';
ae.GridAlpha = .4;
ae.GridLineStyle = '--';
ae.XLabel.String = '[AU]';ae.XLabel.Position(1:2)= [23 0.3];
ae.XAxis.FontSize = labelSize-4;
set(gca,'ytick',[])


set(get(ae,'children'),'clipping','off')% turn off clippings

title('Aliasing index','color',clr.title,'FontSize',labelSize)

%% Save fig
figname = 'fig6_aliasingIndex';
saveas(fig41,[savingDir filesep figname '.png'])
savefig(fig41,[savingDir filesep figname])

%%
fig39 = figure(39);clf
fig39.Position = [423 467 1099 450];
fig39.InvertHardcopy = 'off';
fig39.Color = clr.bck;
ialgo = 1;
hold on

% TF2 = (fft2(MatOutNoInterp{ialgo}));
TF2 = (fft2(MatOut{ialgo}));

TF2_x = smooth(sum(abs(TF2).^2,1)/sum(sum(abs(TF2).^2,1)),1);
TF2_z = smooth(sum(abs(TF2).^2,2)/sum(sum(abs(TF2).^2,1)),1);

peak_z = max(TF2_z(Freq_i_z),[],2);
baseline_z = mean(TF2_z(Freq_base_z),2);
P_alias_z(ialgo) = 20*log10(sum(peak_z./baseline_z));

peak_x = max(TF2_x(Freq_i_x),[],2);
baseline_x = mean(TF2_x(Freq_base_x),2);
P_alias_x(ialgo) = 20*log10(sum(peak_x./baseline_x));

plot(20*log10(TF2_x(1:400)),clr.title,'LineWidth',.2)
xlim([0 325])

for ii=1:2
plot(Freq_i_x(ii,2) + [2 30],20*log10(peak_x(ii))*ones(1,2),':','Color',clr.ax,'LineWidth',1.5)
plot(Freq_i_x(ii,2) + [2 30],20*log10(baseline_x(ii))*ones(1,2),':','Color',clr.ax,'LineWidth',1.5)

text(Freq_i_x(ii,2)+33,20*log10(peak_x(ii))+2,'peak','Color',clr.title,'FontSize',labelSize)
text(Freq_i_x(ii,2)+33,20*log10(baseline_x(ii))-1,'baseline','Color',clr.title,'FontSize',labelSize)

plot(Freq_i_x(ii,2) + [1 1]*15,[20*log10(peak_x(ii)) 20*log10(baseline_x(ii))],'r-','LineWidth',3)

end

grid on
aal = gca;
aal.Color = clr.bck;
aal.XLabel.String = 'Frequency (1/pix)';
aal.YLabel.String = 'Power Specral Density [log]';
aal.XLabel.FontSize = labelSize;
aal.FontSize = labelSize-2;
aal.XLabel.Color = clr.ax;
aal.XAxis.Color = clr.ax;aal.YAxis.Color= clr.ax;
aal.TickDir = 'both';
aal.Box = 'off';
aal.GridAlpha = .4;
aal.GridLineStyle = '--';

%% Save fig
saveas(fig39,[savingDir filesep filename '_AliasingMethod.png'])
savefig(fig39,[savingDir filesep filename '_AliasingMethod'])


%% %% FIG 6 : Saturation curves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure(78);
fig.Position = [841 526 673 359];

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
legend(listAlgoName,'Location','best','FontSize',12)
xlabel('Number of frames')
ylim([0 50])
xticks([0:25:200])
xtickformat('%,.0fk')
ytickformat('%g%%')

grid on

a = gca;
a.Position = [.07 .12 .9 .8];

title('% of nnz pixel')

%%
figname = 'fig6_saturation';
saveas(gcf,[savingDir filesep figname '.png'])
savefig(gcf,[savingDir filesep figname])
% legend('off')
print([savingDir filesep figname 'c'],'-r750','-dpng')

%% %% FIG 2 : Processing Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Build figure
clear ae
fig45 = figure(45);clf
fig45.Position = [530 350 670 260];
fig45.InvertHardcopy = 'off';

ae = axes('Position',[.01 .11 .7 .8]);

TotalTime = sum(ProcessingTime,1)/60;
TotalTime = TotalTime./min(TotalTime);

b_ap = barh(TotalTime,'facecolor',[1 1 1]*0.5,'BarWidth',1);
b_ap.YData(end-1)=7.5;

b_ap.FaceColor = 'flat';
b_ap.CData = ListColor;

text(b_ap.YData-.1,b_ap.XData,string(round(TotalTime,2)),...
    'HorizontalAlignment','right','VerticalAlignment','middle','color',clr.title,'FontSize',labelSize-2)

text(b_ap.YData+.2,b_ap.XData,listAlgoName,...
    'HorizontalAlignment','left','VerticalAlignment','middle','color',clr.title,'FontSize',labelSize)

grid on
ylim([0 8] + [1 -1]*.5)

ae.XAxisLocation = 'origin'; ae.YAxisLocation = 'origin';
ae.XAxis.Color = clr.ax;
ae.YAxis.Color = ae.XAxis.Color;
ae.FontSize = labelSize-2;
ae.XAxis.FontSize = labelSize;
ae.TickDir = 'both';
ae.TickLength = [.005 .0];
ae.Box= 'off';
ae.XLabel.Position(1:2)= [23 0.3];
ae.XAxis.FontSize = labelSize-4;
set(gca,'ytick',[])
ae.XTick = [1:7];

title('Processing time','color',clr.title,'FontSize',labelSize)

%%
figname = 'fig2_processingTiem';
saveas(gcf,[savingDir filesep figname '.png'])
savefig(gcf,[savingDir filesep figname])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Intenstity MatOutnointerp

ibase = 6;
Nloc = cellfun(@(x) sum(x(:)),MatOutNoInterp);

TotalTime = sum(ProcessingTime,1);

%% RADAR DATA
Radar_vivo.Saturation = SaturationEnd;
Radar_vivo.Aliasing = P_alias;
Radar_vivo.Time = TotalTime;
Radar.listName = listAlgoName;
save([workingdir filesep filename '_radar_vivo'],'Radar_vivo')

