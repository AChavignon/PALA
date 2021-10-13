%% PALA_GlobalScore.m : Global Score script for Global Index evaluation
% This script loads results for the 3 datasets (`In Silico PSF`, `In Silico Flow`, `In Vivo Brain`),
% and creates a global score
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
%%
workingdir_flow = [PALA_data_folder '\PALA_data_InSilicoFlow'];
load([workingdir_flow '\PALA_InSilicoFlow_scores.mat'])

workingdir_vivo = [PALA_data_folder '\PALA_data_InVivoRatBrain'];
load([workingdir_vivo '\PALA_InVivoRatBrain_scores.mat'])

workingdir_psf = [PALA_data_folder '\PALA_data_InSilicoPSF'];
Mesh = load([workingdir_psf '\PALA_InSilicoPSF_scores.mat']);

labelSize = 12;
workingdir_results = [PALA_data_folder '\PALA_GlobalScore\'];mkdir(workingdir_results)

%% Create scores for the results of each metric
RMSE = squeeze(Mesh.meanRes(3,:,3)); % at 30 dB
Score.RMSE = 1-RMSE/.5; % [1/1 at 0; 0/11 at .5wv]

Score.Jaccard = Radar.Jaccard/100;
Score.precision = Radar.precision/100;

Gap = Radar.Gap(1:7)/10;
Score.Gap = 1+(-Gap)/.1; % [1/1 at 0; 0/1 at 100um]

ProTime = (Radar_vivo.Time);
Score.Time = 1- log10(ProTime/min(ProTime))/2;clear ProTime

Score.GridIndex = 1-Radar_vivo.GridIndex/30; % [10/10 at 0; 0/10 at -30]

Score.Saturation = Radar_vivo.Saturation; % in percentage

%% Display global score
[ListAlgoName,ListColor,ListMarker,ListShortName] = PALA_GetFormat;

ListScore = [Score.RMSE(:),Score.Jaccard(:),Score.precision(:),Score.Gap(:),Score.Time(:),Score.GridIndex(:),Score.Saturation(:)];
lgdname = {'RMSE','Jaccard','Precision','Gap','Time','Gridding','Saturation'};
lgdshort = {'RMSE','Jacc','P','Gap','Time','Grid','Sat'};
Weight_coef = [2 2 1 1 2 1 1]; % weighting coefficients

% Sort for displaying
PermutValue = [1 5 4 6 2 3 7];
Weight_coef = Weight_coef(PermutValue);
ListScore = ListScore(:,PermutValue);
lgdshort = lgdshort(PermutValue);
lgdname = lgdname(PermutValue);

ListScore = ListScore.*Weight_coef*100./sum(Weight_coef);
GlobalScore = sum(ListScore,2);

fig02 = figure(2);clf;fig02.Position = [450 330 1080 346];axes('Position',[.01 .12 .98 .82])
bb = barh(ListScore,1,'stacked');
for ii=1:numel(bb)
    bb(ii).FaceColor='flat';bb(ii).CData = ListColor(:,:)*1;bb(ii).FaceAlpha = 1-.7*((ii-1)/(numel(bb)-1));
end
ylim([.5 9]);xlim([0 max(GlobalScore)+5])
grid on;aa=gca;set(gca,'ytick',[]);aa.Box= 'off';
text(GlobalScore+1,[1:7], num2str(round(GlobalScore,1)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize)

lgd = legend(bb,lgdname,'Location','northwest','Orientation','horizontal');
txtpos = [0 cumsum(ListScore(end,1:end-1))] + ListScore(end,:)/2;

for ii=1:numel(lgd.String);lgd.String{ii} = [lgd.String{ii} ' (coef ' num2str(Weight_coef(ii)) ')'];end
title('Global scoring')
%%
figname = 'fig2_GlobalScore';
saveas(fig02,[workingdir_results filesep figname '.png'])
print(fig02,[workingdir_results filesep figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% table
Data = [];coltitle = {};
Data(end+1,:) = RMSE;          coltitle{end+1}='RMSE (lambda)';
Data(end+1,:) = Score.RMSE(:)*10;       coltitle{end+1}='Score RMSE (1)';
Data(end+1,:) = Radar_vivo.Time;        coltitle{end+1}='Processing Time (AU)';
Data(end+1,:) = Score.Time(:)*10;       coltitle{end+1}='Score time (1)';
Data(end+1,:) = Gap;                    coltitle{end+1}='Gap (lambda)';
Data(end+1,:) = Score.Gap*10;           coltitle{end+1}='Score gap';
Data(end+1,:) = Radar_vivo.GridIndex;   coltitle{end+1}='Gridding (AU)';
Data(end+1,:) = Score.GridIndex*10;     coltitle{end+1}='Score Gridding (1)';
Data(end+1,:) = Score.Jaccard(:)*10;    coltitle{end+1}='Jaccard (1)';
Data(end+1,:) = Score.precision*10;     coltitle{end+1}='Precision  (1)';
Data(end+1,:) = Radar_vivo.Saturation;  coltitle{end+1}='Saturation (%)';
Data(end+1,:) = Score.Saturation*10;    coltitle{end+1}='Score saturation (1)';
Data(end+1,:) = GlobalScore*10;         coltitle{end+1}='Global (1)';

f = figure(1);clf
uit = uitable(f,'Data',Data,'Units','normalized','Position',[.05 .05 .9 .9]);
uit.ColumnName = ListAlgoName;uit.RowName = coltitle;

figname = 'GlobalScore_values';
save([workingdir_results filesep figname],'Data','coltitle','ListAlgoName','Weight_coef')
fprintf('PALA_GlobalScore.m completed.\n');
