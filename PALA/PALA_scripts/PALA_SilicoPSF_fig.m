%% PALA_SilicoPSF_fig.m : displaying - errors and pairings algorithms for IN SILICO PSF
% Analyses errors of localization for IN SILICO MESH PSF
% Localization positions are loaded and compare to the ground truth. Results are presented
% in different ways to analyses different aspects and behaviors
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
fprintf('Running PALA_SilicoPSF_fig.m - Images generation\n')
workingdir = [PALA_data_folder '\PALA_data_InSilicoPSF'];cd(workingdir)
filename = 'PALA_InSilicoPSF';
myfilepath = [workingdir filesep filename];

load([myfilepath '_sequence.mat'],'Media','filetitle','ll','mesh_x','mesh_z')
myfilepath_fig = [workingdir filesep 'img2' filesep];
if exist(myfilepath_fig)~=7,mkdir(myfilepath_fig);end

[ListAlgoName,ListColor,ListMarker,ShortName] = PALA_GetFormat;
ClutterList = [-60 -40 -30 -25 -20 -15 -10];
load([myfilepath '_LocalMesh'  num2str(30) 'dB']);Nalgo = numel(listAlgo);
clear MatLocFull StaticError

%% Load Data and calculate errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MatPosSim = permute(Media.ListPos(:,[3 1]),[2 3 1]); % true position's of scatterers

ErrorFullAll = [];
for ii=1:numel(ClutterList)
    %% For each dB, load data
    Data = load([myfilepath '_LocalMesh'  num2str(abs(ClutterList(ii))) 'dB'],'MatLocFull','StaticError');disp('Done')
    if abs(ClutterList(ii))==60
        % load the specific localization at -60dB
        StaticError = Data.StaticError(1:2,:);
    end

    % For each frame, calculate the error of localization.
    ErrorFull_i = Data.MatLocFull(:,:,1:size(Media.ListPos,1),:); % remove out of list positions
    ErrorFull_i = ErrorFull_i-MatPosSim;
    % Data are reshape to get a user-friendly organization.
    % for each of the Npoints [21x21 grid] position, we have a 3 random shift sub position (ie [ix/21 + rand/21, iz/21 + rand/21).
    % The calcultation is computed fo each Algo (actually 7), and for Nt(2)  occasion of noise.
    % --> for each alog, we have 21x21x3x2 = 2646 frames

    ErrorFull_i = reshape(permute(ErrorFull_i,[1,3,4,2]),2,[],Nalgo);
    ErrorFull_i = ErrorFull_i - permute(StaticError(1:2,:),[1 3 2]);
    ErrorFull_i(3,:,:) = sqrt(sum(ErrorFull_i(1,:,:).^2+ErrorFull_i(2,:,:).^2,1));
    % Final data with all errors values:
    %       Dim 1 : localization error [z x]
    %       Dim 2 : frame 21x21 x 3 x 2
    %       Dim 3 : list algo (7)
    %       Dim 4 : list of dB
    ErrorFullAll(:,:,:,ii)=ErrorFull_i;
end
clear Data ErrorFull_i

%% FIG 4 : diagrams Errors Mesh lateral/axial/RMSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanRes = squeeze(mean(ErrorFullAll,2));varRes = squeeze(std(ErrorFullAll,0,2));

fig45 = figure(45);clf;fig45.Position = [297 189 800 535];
aa = tight_subplot(3,1,.07,0.02,[.05 .005]);
ii=1;axes(aa(ii));
PALA_boxplot_multi([],ErrorFullAll(2,:,:,:),12,ListColor',9,0);
title('Axial Error','Position',[4 .4 0])

for ii=1:numel(ClutterList)
    tt = text(.55+(ii-1)/(numel(ClutterList)-1)*.9,-.6,[num2str(abs(ClutterList(ii)))],'VerticalAlignment','bot','HorizontalAlignment','center','FontSize',7);
end
tt.String = [tt.String 'dB'];
grid on;a = gca;a.FontSize = 10;
a.YTick = [-.5:.25:.5];a.YTickLabelMode = 'auto';
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;ylim([-.5 .5])

ii=2;axes(aa(ii));
PALA_boxplot_multi([],ErrorFullAll(1,:,:,:),15,ListColor',9,0);
title('Lateral Error','Position',[4 .4 0])
grid on;a = gca;a.FontSize = 10;
a.YTick = [-.5:.25:.5];a.YTickLabelMode = 'auto';
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;ylim([-.5 .5])

ii=3;axes(aa(ii));
PALA_boxplot_multi([],ErrorFullAll(3,:,:,:),15,ListColor',9,0);
title('RMSE','Position',[4 .5 0])
grid on;a = gca;a.FontSize = 10;
a.YTick = [0:.25:1];a.YTickLabelMode = 'auto';
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;ylim([0 .75])

aa(3).XTickLabel = ListAlgoName;
save([myfilepath '_scores.mat'],'NoiseParam','listAlgo','ULM','ClutterList','meanRes')

%% Save fig
figname = 'fig4_ErrorsMesh';
print(fig45,[myfilepath_fig figname],'-dpng','-r300')
print(fig45,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% FIG 4 : Error maps at 30dB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The value in each pixel represents the absolute mean error of localization for bubbles in that pixel.
fig46 = figure(46);fig46.Position = [510 521 600 466];clf;colormap jet

Caxis.x = [0 1]*.3;Caxis.y = [0 1]*.3;Caxis.RMSE = [0 1]*.4;
Algo2Disp = [2 5 6 7]; % Display only few algo for paper
aa = tight_subplot(3,numel(Algo2Disp),[.03 .01],[.02 .05],[.07 .07]);

ErrorMaps = ErrorFullAll(:,:,:,find(ClutterList==-30));
ErrorMaps = reshape(ErrorMaps,3,numel(ll),numel(ll),[],Nalgo);
ErrorMaps = squeeze(mean(ErrorMaps,4));
ErrorMaps = permute(ErrorMaps,[2 3 1 4]);
ErrorMaps = abs(ErrorMaps); % [Nx Nz Npos Nalgo]
% For each algorithm, display the error map corresponding to the 21x21 grid in z, x, and RMSE.
% The value is averaged with the number of occurrence of noise and the number of random shift (3x2 per pixel)

for ii = 1:numel(Algo2Disp)
    ialgo = Algo2Disp(ii);
    axes(aa(ii))
    imagesc(ll,ll,ErrorMaps(:,:,2,ialgo)),axis image
    caxis(Caxis.x)
    title(ListAlgoName(ialgo))
    if ii==1
        ylabel('Lateral Error');axis on;
    else
        set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
end
clb = colorbar;
clb.Position(1) = aa(numel(Algo2Disp)).Position(1)+aa(numel(Algo2Disp)).Position(3)*1.05;clb.Position([2 4]) = aa(numel(Algo2Disp)).Position([2 4]);

for ii = 1:numel(Algo2Disp)
    ialgo = Algo2Disp(ii);
    axes(aa(ii+numel(Algo2Disp)))
    imagesc(ll,ll,ErrorMaps(:,:,1,ialgo)),axis image
    caxis(Caxis.y)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    if ii==1,ylabel('Axial Error');end
end
cbr_z = colorbar;
cbr_z.Position(1)=aa(numel(Algo2Disp)*2).Position(1)+aa(numel(Algo2Disp)*2).Position(3)*1.05;cbr_z.Position([2 4])=aa(numel(Algo2Disp)*2).Position([2 4]);

for ii = 1:numel(Algo2Disp)
    ialgo = Algo2Disp(ii);
    axes(aa(ii+numel(Algo2Disp)*2))
    imagesc(ll,ll,ErrorMaps(:,:,3,ialgo)),axis image
    caxis(Caxis.RMSE)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    if ii==1,ylabel('RMSE');end
end
cbr_r = colorbar;
cbr_r.Position(1)=aa(numel(Algo2Disp)*3).Position(1)+aa(numel(Algo2Disp)*3).Position(3)*1.05;cbr_r.Position([2 4])=aa(numel(Algo2Disp)*3).Position([2 4]);

%% Save fig
figname = 'fig4_ErrorsMaps_30dB';
print(fig46,[myfilepath_fig figname],'-dpng','-r300')
print(fig46,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% FIG 4 : Display profil %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = load([myfilepath '_IQ001.mat']);
IQ = abs(temp.IQ(:,:,1:size(temp.Media.ListPos,1)));
IQ = IQ./max(IQ(:));

listPost = temp.ListPos(1:size(IQ,3),:);
listPost2 = listPost(:,1:3)-temp.Media.CenterPoint;
listPost2 = listPost2./temp.Media.Displacement;
listPost2 = round(listPost2);
listPost2 = listPost2(1:21*21,:);

rdisp = [-2:2];
listX = [-9:3:9];nx = numel(listX);
listY = [-9:3:9];ny = numel(listY);

m_w = .02;m_h = .09;gap_w = .008;gap_h = .1;
pr_s_w = (1-2*m_w-gap_w*(nx*2-1))/nx/2;pr_s_hs = .15;
pr_s_h = (1-m_h*2-gap_h-pr_s_hs)/1;

posax = @(ix,ip,isc)[m_w+(ix-1)*pr_s_w+(ix-1)*gap_w,...
    m_h+ip*pr_s_h+isc*pr_s_hs+(ip+isc)*gap_h,...
    pr_s_w,...
    pr_s_hs*(mod(ip+isc,3)==1) + pr_s_h*(1-(mod(ip+isc,3)==1))];

half_fontsize = 9;
f1 = figure(1);clf;f1.Position = [46 600 1479 150];
for ix = 1:nx
    ind =and(listPost2(:,1)==listX(ix),listPost2(:,3)==0);ind = find(ind);
    axes('Position',posax(ix,1,0))
    plot(listX(ix)*Media.Displacement(1),0,'k.','MarkerSize',10);axis image
    xlim([-1 1]*.5);ylim([-1 1]*.5)
    aa=gca;aa.XAxisLocation = 'origin';aa.YAxisLocation = 'origin';aa.XAxis.Color = [1 1 1]*.5;aa.YAxis.Color = aa.XAxis.Color;
    if ix==1
        text(-.5,0,'-.5','FontSize',half_fontsize,'VerticalAlignment','middle','HorizontalAlignment','right');
        text(.5,0,'.5','FontSize',half_fontsize,'VerticalAlignment','middle','HorizontalAlignment','left');
        text(0,-.5,'-.5','FontSize',half_fontsize,'VerticalAlignment','top','HorizontalAlignment','center');
        text(0,.5,'.5','FontSize',half_fontsize,'VerticalAlignment','bot','HorizontalAlignment','center');
    end
    % lateral profil
    ap = axes('Position',posax(ix,0,0));hold on
    ap.XAxisLocation = 'origin';ap.YAxisLocation = 'origin';
    p0 = plot(rdisp,IQ(rdisp+5,6,ind(1)),'x-','LineWidth',1.2,'MarkerSize',6,'DisplayName','Axial profile');
    p1 = plot(rdisp,IQ(5,rdisp+6,ind(1)),'LineWidth',p0.LineWidth,'MarkerSize',p0.MarkerSize,'Marker','+','DisplayName','Lateral profile');
    plot(listX(ix)*Media.Displacement(1)*[1 1],[0 1],'--','Color',p1.Color,'LineWidth',1.5,'DisplayName','Real position');
    grid on;ylim([0 1])
    set(gca,'YTickLabel',{})
    if ix>1;set(gca,'XTickLabel',{});end
    ap.XAxis.TickLength=[.08,0];ap.XAxis.MinorTick = 'on';ap.XAxis.MinorTickValues=[-2:.5:2];
end

for ix = 1:nx
    ind = and(listPost2(:,1)==0,listPost2(:,3)==listY(ix));ind = find(ind);
    axes('Position',posax(ix+nx,1,0))
    plot(0,listY(ix)*Media.Displacement(1),'k.','MarkerSize',10)
    axis image
    xlim([-1 1]*.5);ylim([-1 1]*.5)
    aa=gca;aa.XAxisLocation = 'origin';aa.YAxisLocation = 'origin';
    aa.XAxis.Color = [1 1 1]*.5;aa.YAxis.Color = aa.XAxis.Color;
    if ix==1
        text(-.5,0,'-.5','FontSize',half_fontsize,'VerticalAlignment','middle','HorizontalAlignment','right');
        text(.5,0,'.5','FontSize',half_fontsize,'VerticalAlignment','middle','HorizontalAlignment','left');
        text(0,-.5,'-.5','FontSize',half_fontsize,'VerticalAlignment','top','HorizontalAlignment','center');
        text(0,.5,'.5','FontSize',half_fontsize,'VerticalAlignment','bot','HorizontalAlignment','center');
    end
    % Axial profil
    ap=axes('Position',posax(ix+nx,0,0));hold on
    ap.XAxisLocation = 'origin';ap.YAxisLocation = 'origin';
    p0 = plot(rdisp,IQ(rdisp+5,6,ind(1)),'x-','LineWidth',1.2,'MarkerSize',6,'DisplayName','Axial profile');
    p1=plot(rdisp,IQ(5,rdisp+6,ind(1)),'LineWidth',p0.LineWidth,'MarkerSize',p0.MarkerSize,'Marker','+','DisplayName','Lateral profile');
    plot(listY(ix)*Media.Displacement(1)*[1 1],[0 1],'--','Color',p0.Color,'LineWidth',1.5,'DisplayName','Real position');

    grid on;ylim([0 1])
    if ix>0;set(gca,'XTickLabel',{});set(gca,'YTickLabel',{});end
    ap.XAxis.TickLength=[.08,0];ap.XAxis.MinorTick = 'on';ap.XAxis.MinorTickValues=[-2:.5:2];
end

%% Save fig
figname = 'fig4_ProfilsIQ3_c';
print(f1,[myfilepath_fig figname '.eps'],'-depsc2','-tiff','-r300','-painters')
print(f1,[myfilepath_fig figname],'-dpng','-r300')

fprintf('PALA_SilicoPSF_fig.m done.\n');
