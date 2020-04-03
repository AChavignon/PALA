%% ULM PARLA: Global score
% This script loads results for the 3 datasets (InSilicoPSF, InSilicoFlow, InVivoBrain),
% and creates a global score
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
%%
workingdir = [PARLA_data_folder '\PARLA_data_InSilicoFlow'];
% filename = '2D_flow_1203-180447';
load([workingdir '\2D_flow_1203-180447_radar_simu.mat'])

workingdir_vivo = [PARLA_data_folder '\PARLA_data_InVivoBrain\PARLA_InVivoRatBrain_'];
load([workingdir_vivo 'radar_vivo'])

workingdir_psf = [PARLA_data_folder '\PARLA_data_InSilicoPSF'];
Mesh = load([workingdir_psf '\2D_mesh_0205_120154Mesh_scoring.mat'],'meanRes','listAlgoName2','listAlgo','listdB','NoiseParam');

labelSize = 12;
myfilepath_res = [workingdir filesep 'img_paper4' filesep];

%% Create socre for results
RMSE = squeeze(Mesh.meanRes(3,:,3));
RMSE_score = 1+(0-RMSE)/.5; % [100 at 0; 0 at .5pix]

figure(1);clf
subplot 211;plot(RMSE);subplot 212;plot(RMSE_score)
Score.RMSE = RMSE_score;

Score.Jaccard = Radar.Jaccard/100;
Score.precision = Radar.precision/100;

Gap = Radar.Gap(1:7)/10;
Gap_score = 1+(0-Gap)/1; % [100 at 0; 0 at 1pix]

figure(1);clf
subplot 211;plot(Gap);subplot 212;plot(RMSE_score)
Score.Gap = Gap_score;


ProTime = (Radar_vivo.Time);%ProTime = ProTime./min(ProTime);
ProTime = 20*log10(ProTime);
Score.Time = 1-(min(ProTime)-ProTime)/-40;


Aliasing = Radar_vivo.Aliasing;
Score.Aliasing = 1 - (0-Aliasing)/-30; % [100 at 0; 0 at -30]

Saturation = Radar_vivo.Saturation;
Score.Saturation = Saturation; % in percentage

%% PLOT SCORRING BAR
[ListAlgoName,ListColor,ListMarker] = GetFormatFightClub;

PermutValue = [1 5 4 6 2 3 7];

y = [Score.RMSE(:),Score.Jaccard(:),Score.precision(:),Score.Gap(:),Score.Time(:),Score.Aliasing(:),Score.Saturation(:)];
lgdname = {'RMSE','Jaccard','Precision','Gap','Time','Aliasing','Saturation'};
lgdshort = {'RMSE','Jacc','P','Gap','Time','Alias','Sat'};

coef = [2 2 1 1 2 1 1];
y = y.*coef;

coef = coef(PermutValue);
y = y(:,PermutValue);
lgdshort = lgdshort(PermutValue);
lgdname = lgdname(PermutValue);

GlobalScore = sum(y,2);

fig02 = figure(2);clf;
fig02.Position = [450 330 1080 346];
axes('Position',[.01 .12 .98 .82])

bb = barh(y*10,1,'stacked');
for ii=1:numel(bb)
    bb(ii).FaceColor='flat';
    bb(ii).CData = ListColor(:,:)*1;
    bb(ii).FaceAlpha = 1-.7*((ii-1)/(numel(bb)-1));
end

xlabel('Scoring')
ylim([.5 9])
xlim([0 max(GlobalScore*10)+5])
grid on
aa=gca;
set(gca,'ytick',[])
aa.Box= 'off';

text(GlobalScore*10+1,[1:7],ListAlgoName,...
    'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize)

lgd =legend(bb,lgdname,'Location','northwest','Orientation','horizontal');
txtpos = [0 cumsum(y(end,1:end-1))] + y(end,:)/2;
% text(txtpos*10,7*ones(1,numel(txtpos)),lgdshort,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8)

for ii=1:numel(lgd.String)
    lgd.String{ii} = [lgd.String{ii} ' (coef ' num2str(coef(ii)) ')'];
end

title('Global scoring')
%%
figname = 'fig2_GlobalScore';
saveas(fig02,[myfilepath_res filesep figname '.png'])
savefig(fig02,[myfilepath_res filesep figname])



