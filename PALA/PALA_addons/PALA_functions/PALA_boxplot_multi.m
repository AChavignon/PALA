function val = PALA_boxplot_multi(label_in,val_Err,labelSize,InPutColor,txtsize,dispVal)
%% function val = PALA_boxplot_multi(label_in,val_Err,labelSize,InPutColor,txtsize,dispVal)
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

if isempty(labelSize),labelSize = 12;end
% Compute statistics

if iscell(val_Err)
    ival = size(val_Err,2);
    for i1=1:size(val_Err,1)
        for i2=1:size(val_Err,2)
            val_med(i1,i2) = median(val_Err{i1,i2});
            val_mean(i1,i2) = mean(val_Err{i1,i2});
            val_std(i1,i2) = std(val_Err{i1,i2});
            val_d1(i1,i2) = quantile(val_Err{i1,i2},.05);
            val_d9(i1,i2) = quantile(val_Err{i1,i2},.95);
            val_q1(i1,i2) = quantile(val_Err{i1,i2},.25);
            val_q3(i1,i2) = quantile(val_Err{i1,i2},.75);
        end
    end
else
    val_Err = squeeze(val_Err);

    % Statistics
    ival = size(val_Err,3);
    val_mean = mean(val_Err,1);
    val_std = std(val_Err,0,1);
    val_med = median(val_Err,1);
    val_d1 = quantile(val_Err,.05,1);
    val_d9 = quantile(val_Err,.95,1);
    val_q1 = quantile(val_Err,.25,1);
    val_q3 = quantile(val_Err,.75,1);
end

val = cat(3,val_mean,val_std);
%% Displaying
cla
hold on
box_width = .94/ival;
spacing = box_width/2*(ival-1);

posX = linspace(-spacing,spacing,ival);

val_mean = squeeze(val_mean);
val_med = squeeze(val_med);val_std = squeeze(val_std);
val_d1 = squeeze(val_d1);val_d9 = squeeze(val_d9);
val_q1 = squeeze(val_q1);val_q3 = squeeze(val_q3);

for ii=1:size(val_mean,1)
    for ia = 1:size(val_mean,2)
        col = InPutColor(:,ii)*.8;
        col = col * ( 1+0.5*(ia-1)/size(val_mean,2));
        col(col>1)=1;
%         col(4) = 1-0.8*(ia-1)/size(val_mean,2);

        p=plot(ii + [0 0] + posX(ia),[val_d1(ii,ia) val_d9(ii,ia)],'.-k','linewidth',1);
        p.Marker = '.';p.MarkerSize = 10;

        posrect = [ii-box_width/2,val_q1(ii,ia),box_width,val_q3(ii,ia)-val_q1(ii,ia)];
        posrect(1) = posrect(1)+ posX(ia);

        rectangle('position',posrect,'edgecolor','k','FaceColor',col,'linewidth',.5)
        plot(ii + [-1 1].*box_width/2+posX(ia),val_med(ii,ia)*[1 1],'k-','linewidth',1)
        plot(ii+posX(ia),val_mean(ii,ia),'ko','MarkerFaceColor','k','MarkerSize',4)

%         posrect = [ii-box_width/2,val_mean(ii,ia)-val_std(ii,ia),box_width/2,val_std(ii,ia)*2];
%         posrect(1) = posrect(1)+posX(ia);
%         rectangle('position',posrect,'edgecolor','none','FaceColor',col*.9)

        if dispVal
        text(ii+box_width/2 +posX(ia)-.1,val_d9(ii,ia)+.02,...
            [num2str(round(val_mean(ii,ia),2)) '\pm' num2str(round(val_std(ii,ia),2),'%.2f')],...
            'fontsize',txtsize,'interpreter','tex','HorizontalAlignment','left','VerticalAlignment','middle','color','k','Rotation',90)
        end
    end

    xlim([0.5 size(val_mean,1)+.5]);
    xticks([1:size(val_mean,1)]);    xticklabels({})
    ax = gca;ax.XAxis.Color='k';ax.YAxis.Color='k';ax.FontSize = labelSize;
end

end
