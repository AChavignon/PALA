function PALA_boxplot(label_in,val_Err,varargin)
%% function PALA_boxplot(label_in,val_Err)
% Custom function to display box plot.
%
% Created by Arthur Chavignon 2019
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

labelSize=10;

%% Statistics
if iscell(val_Err)
    val_med = cellfun(@median,val_Err);
    val_mean = cellfun(@mean,val_Err);
    val_std = cellfun(@std,val_Err);
    val_d1 =  cellfun(@(x) quantile(x,.05,1),val_Err);
    val_d9 = cellfun(@(x) quantile(x,.95,1),val_Err);
    val_q1 =  cellfun(@(x) quantile(x,.25,1),val_Err);
    val_q3 =  cellfun(@(x) quantile(x,.75,1),val_Err);
else
    val_med = mean(val_Err);
    val_mean = mean(val_Err);
    val_std = std(val_Err);
    val_d1 = quantile(val_Err,.05,1);
    val_d9 = quantile(val_Err,.95,1);
    val_q1 = quantile(val_Err,.25,1);
    val_q3 = quantile(val_Err,.75,1);
end

tmp = strcmpi(varargin,'InPutColor');
if any(tmp), InPutColor= varargin{find(tmp)+1}; else, InPutColor=ones(numel(val_mean),3)*.5;end
tmp = strcmpi(varargin,'box_width');
if any(tmp), box_width= varargin{find(tmp)+1}; else, box_width=.4;end

%%
hold on

for ii=1:numel(val_mean)
    p=plot(ii + [0 0],[val_d1(ii) val_d9(ii)],'-k','linewidth',1.5);
    p.Marker = '.';p.MarkerSize = labelSize;

    posrect = [ii-box_width/2,val_q1(ii),box_width,val_q3(ii)-val_q1(ii)];
    rectangle('position',posrect,'edgecolor','k','FaceColor',InPutColor(ii,:),'linewidth',1)
    plot(ii + [-1 1].*box_width/2,val_med(ii)*[1 1],'k--','linewidth',1)
    plot(ii,val_mean(ii),'ko','MarkerFaceColor','k','MarkerSize',4)

    text(ii+box_width/2 + .02,double(val_mean(ii)),...
        [num2str(round(val_mean(ii),3)) '\pm' num2str(round(val_std(ii),3),'%.2f')],...
        'fontsize',labelSize,'interpreter','tex','HorizontalAlignment','left','VerticalAlignment','middle','color','k','Rotation',80)
end

xlim([0.5 numel(val_mean)+.5]);xticks([1:numel(val_mean)]);xticklabels({})
ax = gca;ax.XAxis.Color='k';ax.YAxis.Color='k';ax.FontSize = labelSize;
end
