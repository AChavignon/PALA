function myBoxPlotFC(label_in,val_Err)
%% function myBoxPlotFC(label_in,val_Err)
% Custom function to display box plot.
%
% Created by Arthur Chavignon 2019
%
% DATE 2020.03.30 - VERSION 1.0.0
% AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
% Code Available under Creative Commons Non-Commercial 4.0
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
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

%%
hold on
box_width = .3;
  
for ii=1:numel(val_mean)
    p=plot(ii + [0 0],[val_d1(ii) val_d9(ii)],'-k','linewidth',1.5);
    p.Marker = '.';p.MarkerSize = labelSize;

    posrect = [ii-box_width/2,val_q1(ii),box_width,val_q3(ii)-val_q1(ii)];

    rectangle('position',posrect,'edgecolor','k','FaceColor',[1 1 1]*.5,'linewidth',1)
    plot(ii + [-1 1].*box_width/2,val_med(ii)*[1 1],'k--','linewidth',1)

    plot(ii,val_mean(ii),'ko','MarkerFaceColor','k','MarkerSize',4)

    text(ii+box_width/2 + .02,val_mean(ii),...
        [num2str(round(val_mean(ii),3)) '\pm' num2str(round(val_std(ii),3),'%.2f')],...
        'fontsize',labelSize,'HorizontalAlignment','left','VerticalAlignment','middle','color','k')    
end

xlim([0.5 numel(val_mean)+.5]);xticks([1:numel(val_mean)]);xticklabels({})
ax = gca;ax.XAxis.Color='k';ax.YAxis.Color='k';ax.FontSize = labelSize;

end