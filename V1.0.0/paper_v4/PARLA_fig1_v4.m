%% ULM PARLA: materials for figure 1
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
PARLA_data_folder = 'F:\ArthurC\Data\Fight_Club\Data';
addpath(genpath(PARLA_addons_folder))
savinfig1 = [PARLA_data_folder '\Fig1'];
%% IN VIVO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
workingdir = [PARLA_data_folder];
cd(PARLA_data_folder)
%% Matout Vivo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_invivo = [PARLA_data_folder filesep '\PARLA_data_InVivoBrain\'];

temp = load([path_invivo 'fightclub_round2/PARLA_InVivoRatBrain_MatOut_multi0' '.mat']);

figure(1);clf;set(gcf,'Position',[680 585 598 393]);
axes('Position',[0 0 1 1]);
matOut = temp.MatOut{end}(20:end-10,15:end-10).^(1/3);
imagesc(matOut)
colormap hot
caxis([0 10])
axis image
axis off

figname = 'fig1_matout_radial_vivo';
saveas(gcf,[savinfig1 filesep figname '.png'])
savefig(gcf,[savinfig1 filesep figname])
WriteTif(matOut,hot,[savinfig1 filesep figname '.tif'],'caxis',[0 10],'overwrite',1)

clear temp
%% IQ no filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IQ =  fastLoad([workingdir '\IQ\SLFonc7_vib_d_8_bulles_IQ104.mat.fs']);
load([PARLA_data_folder '\IQ\PARLA_InVivoRatBrain_104.mat'],'IQ')
IQ = IQ(2:end-1,2:118,:);
dB = IQdB(IQ);

figure(2);clf;set(gcf,'Position',[680 585 598 393]);
axes('Position',[0 0 1 1]);
imagesc(dB(:,:,12),[-40 0])
colormap gray
axis image
axis off

figname = 'fig1_IQnofilt_at40dB';
saveas(gcf,[savinfig1 filesep figname '.png'])
savefig(gcf,[savinfig1 filesep figname])
WriteTif(dB(:,:,12),gray,[savinfig1 filesep figname '.tif'],'caxis',[-40 0],'overwrite',1)


%% IQ filtered %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQfilt = SVDfilter(IQ,[10 800]);
dBfilt = IQdB(IQfilt);

figure(3);clf;set(gcf,'Position',[680 585 598 393]);
axes('Position',[0 0 1 1]);
imagesc(dBfilt(:,:,12),[-40 0])
colormap gray
axis image
axis off
 
figname = 'fig1_IQfilt_at40dB';
saveas(gcf,[savinfig1 filesep figname '.png'])
savefig(gcf,[savinfig1 filesep figname])
WriteTif(dBfilt(:,:,12),gray,[savinfig1 filesep figname '.tif'],'caxis',[-40 0],'overwrite',1)

%%
clear IQfilt IQ dBfilt dB 
clearvars -except savinfig1

%% MESh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% IQ noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
workingdir = [PARLA_data_folder '\PARLA_data_InSilicoPSF'];
% filename = '2D_mesh_0205_120154';];
cd(workingdir)
temp = load('2D_mesh_0205_120154_IQ001','Media','IQ');
load('2D_mesh_0205_120154_LocalMesh10dB','NoiseParam');
load('2D_mesh_0205_120154Mesh_scoring','listdB');

IQ = abs(temp.IQ(:,:,100:110));
figname = 'fig1_mesh_snr';

figure(4);clf
set(gcf,'Position',[680 818 1030 160])
for ii=1:numel(listdB)
    NoiseParam.clutterdB = listdB(ii);
%     subplot(1,numel(listdB),ii)
    axes('Position',[.01+.95/7*(ii-1),.1,.95/7,.8])
    IQ_speckle = IQ+imgaussfilt(max(IQ(:))*10^(NoiseParam.clutterdB/20) + reshape(wgn(numel(IQ),1,NoiseParam.Power,NoiseParam.Impedance),size(IQ,1),size(IQ,2),[])*max(IQ(:))*10^((NoiseParam.amplCullerdB+NoiseParam.clutterdB)/20),NoiseParam.SigmaGauss);
    IQ_speckle = IQdB(IQ_speckle);
    imagesc(IQ_speckle(:,:,2),[-40 0])
    axis image; axis off; colormap gray
    title([num2str(NoiseParam.clutterdB) 'dB'])
    
    WriteTif(IQ_speckle(:,:,2),gray,[savinfig1 filesep figname num2str(abs(listdB(ii))) 'dB.tif'],'caxis',[-40 0],'overwrite',1)
end
clb = colorbar;clb.Position = [.96 .1 .01 .8];

% saveas(gcf,[savinfig1 filesep figname '.png'])
% savefig(gcf,[savinfig1 filesep figname])

%% Target pin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Media = temp.Media;
figure(5);clf
set(gcf,'Position',[500 500 300 300])
spa = 8;
plot(Media.ListPos(1:spa:441,1)-Media.CenterPoint(1),Media.ListPos(1:spa:441,3)-Media.CenterPoint(3),'ko','MarkerSize',4,'MarkerFaceColor','k')
axis image
xlim([-.5 .5]);ylim([-.5 .5])

grid on
aa=gca;aa.Position = [.05 .05 .9 .9];
aa.XAxisLocation = 'origin';aa.YAxisLocation = 'origin';
aa.XTickLabel = {};aa.YTickLabel = {};
aa.XAxis.Color = [1 1 1]*.6;
aa.YAxis.Color = aa.XAxis.Color;

figname = 'fig1_mesh_target';
saveas(gcf,[savinfig1 filesep figname '.png'])
saveas(gcf,[savinfig1 filesep figname '.svg'])
savefig(gcf,[savinfig1 filesep figname])

%% IN SILICO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Flow SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
workingdir = [PARLA_data_folder '\PARLA_data_InSilicoFlow'];
cd(workingdir)
temp = load([workingdir '\IQ\2D_flow_1203-180447_IQ004.mat']);
load([workingdir '\config2\2D_flow_1203-180447_Stats_multi10dB.mat'],'ClutterList','NoiseParam')
IQ = abs(temp.IQ(:,:,1:110));

figname = 'fig1_flow_snr';

figure(6);clf
set(gcf,'Position',[286 300 1472 581])
aa = tight_subplot(2,3,.01,0,0);
for ii=1:numel(ClutterList)
    NoiseParam.clutterdB = ClutterList(ii);
%     subplot(2,4,ii)
    axes(aa(ii))
    IQ_speckle = IQ+imgaussfilt(max(IQ(:))*10^(NoiseParam.clutterdB/20) + reshape(wgn(numel(IQ),1,NoiseParam.Power,NoiseParam.Impedance),size(IQ,1),size(IQ,2),[])*max(IQ(:))*10^((NoiseParam.amplCullerdB+NoiseParam.clutterdB)/20),NoiseParam.SigmaGauss);
    IQ_speckle = IQdB(IQ_speckle);
    imagesc(IQ_speckle(:,:,20),[-40 0])
    axis image; axis off; colormap gray
    
    WriteTif(IQ_speckle(:,:,20),gray,[savinfig1 filesep figname num2str(abs(ClutterList(ii))) 'dB.tif'],'caxis',caxis,'overwrite',1)
end
%%
saveas(gcf,[savinfig1 filesep figname '.png'])
savefig(gcf,[savinfig1 filesep figname])

%% MatOut Target %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = load([workingdir '\config2\2D_flow_1203-180447_MatOut_multi_30dB' '.mat']);

figure(1);clf;set(gcf,'Position',[680 619 611 360]);
axes('Position',[0 0 1 1]);
matOut = imgaussfilt(temp.MatOut{end}(:,:).^(1),2);
imagesc(matOut)

colormap hot
caxis([0 40])
axis image
axis off
%%
figname = 'fig1_matout_silico';
% saveas(gcf,[savinfig1 filesep figname '.png'])
savefig(gcf,[savinfig1 filesep figname])
WriteTif(matOut,hot,[savinfig1 filesep figname '.tif'],'caxis',caxis)
% clear temp

%% Display streamline
cd(workingdir)
load([workingdir '\2Dpaper3_config.mat']);
load([workingdir '\2Dpaper3_positions_MediaPos']);

raio = 1.6;
figure(7);clf
pos = [557 230 1012 675];
pos(3:4) = pos(3:4).*[0.7750 0.8150];
pos = [10 40 784.3000 465.1250];pos(3:4) = pos(3:4)*raio;
set(gcf,'Position',pos)
set(gcf,'Color','w')
axes('Position',[0 0 1 1])
hold on,axis image

for iCurve = 1:numel(H)
    pp=fnplt(H{iCurve}.curve);
    rr = H{iCurve}.d;
    rr = max(rr,.2);
    plot(pp(1,:),pp(3,:),'LineWidth',rr*4*raio,'Color',[1 1 1]*.7)
end

CurAx = gca;
CurAx.YDir = 'reverse';
% CurAx.View = [0 0];
xlim(PData(1).Origin(1) + [0 PData(1).Size(2)*PData(1).PDelta(1)] )
ylim([P.startDepth P.endDepth])

for ii=[800 1200]
mm = MediaPos{ii};
plot(mm(:,1),mm(:,3),'ko','MarkerSize',4,'MarkerFaceColor','k')
end
axis off
%%
figname = 'fig1_target_silico';
% saveas(gcf,[savinfig1 filesep figname '.tif'])
saveas(gcf,[savinfig1 filesep figname '.svg'])
savefig(gcf,[savinfig1 filesep figname])




