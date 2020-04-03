%% ULM PARLA: displaying - errors and pairings algorithms for IN SILICO PSF
% Analys errors of localization for IN SILICO MESH PSF
% Localization positions are loaded and compare to the ground truth. Results are presented
% in different ways to analys different aspects and behaviours
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
%% Select IQ and Media folder
workingdir = [PARLA_data_folder '\PARLA_data_InSilicoPSF'];
filename = '2D_mesh_0205_120154';
cd(workingdir)

myfilepath = [workingdir filesep filename];
mkdir([workingdir filesep 'img'])
resultFile = [workingdir filesep 'img' filesep filename];

load([myfilepath '_sequence.mat'],'Media','filetitle','ll','mesh_x','mesh_z')
suffixe = '1';

% savingDir = 'D:\ArthurC\Fight_Club\PARLA_figures';
myfilepath_res = [workingdir filesep 'img_paper4' filesep];
if exist(myfilepath_res)~=7
    mkdir(myfilepath_res)
end

[ListAlgoName,ListColor,ListMarker,ShortName] = GetFormatFightClub;
ClutterList = [-60 -40 -30 -25 -20 -15 -10];
load([myfilepath '_LocalMesh'  num2str(30) 'dB']);
clear MatLocFull StaticError
Nalgo = numel(listAlgo);

%% Load Data and calculate errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MatPosSim = permute(Media.ListPos(:,[3 1]),[2 3 1]); % theorical position's of scatterers

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


%% %% FIG 4 : diag Errors Mesh lateral Axial RMSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lblsize = 12;
meanRes = squeeze(mean(ErrorFullAll,2));
varRes = squeeze(std(ErrorFullAll,0,2));

fig45 = figure(45);clf;fig45.Position = [297 189 800 535];

aa = tight_subplot(3,1,.07,0.02,[.05 .005]);
ii=1;
axes(aa(ii));
myBoxPlotFC_multi([],ErrorFullAll(ii,:,:,:),12,ListColor',9,0);
% ylabel('Axial error')

for ii=1:numel(ClutterList)
    tt = text(.55+(ii-1)/(numel(ClutterList)-1)*.9,-.6,[num2str(ClutterList(ii))],'VerticalAlignment','bot','HorizontalAlignment','center','FontSize',7);
end
tt=text(.55+(ii)/(numel(ClutterList)-1)*.9,-.6,'dB','VerticalAlignment','bot','HorizontalAlignment','center','FontSize',7);
grid on
a = gca;a.FontSize = 10;
a.YTick = [-.5:.25:.5];a.YTickLabelMode = 'auto';
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;
ylim([-.5 .5])

ii=2;
axes(aa(ii));
myBoxPlotFC_multi([],ErrorFullAll(ii,:,:,:),15,ListColor',9,0);
% ylabel('Lateral error')
grid on;a = gca;a.FontSize = 10;
a.YTick = [-.5:.25:.5];a.YTickLabelMode = 'auto';
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;
ylim([-.5 .5])

ii=3;
axes(aa(ii));
myBoxPlotFC_multi([],ErrorFullAll(ii,:,:,:),15,ListColor',9,0);
% ylabel('RMSE')
grid on;a = gca;a.FontSize = 10;
a.YTick = [0:.25:1];a.YTickLabelMode = 'auto';
a.YMinorGrid = 'on';a.MinorGridAlpha = .15;
ylim([0 .75])

% aa(3).XTickLabel = ListAlgoName;
% suptitle(['Localization errors for [' num2str(listdB) '] dB (' num2str(size(ErrorFull,3)*size(ErrorFull,4)) ' points per configuration)'])

%% Save fig
figname = 'fig4_ErrorsMesh';
savefig(fig45,[myfilepath_res filesep figname])
print(fig45,[myfilepath_res filesep figname],'-dpng','-r600')

%% %% FIG 4 : Error maps 30dB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The value in each pixel represents the absolute mean error of localization for bubbles
% in that pixel.

fig46 = figure(46);fig46.Position = [510 521 585 466];clf
colormap jet

Caxis.x = [0 1]*.2;
Caxis.y = [0 1]*.3;
Caxis.RMSE = [0 1]*.4;

Algo2Disp = [2 5 6 7]; % Display only few algo
aa = tight_subplot(3,numel(Algo2Disp),[.03 .01],[.02 .05],[.04 .07]);

ErrorMaps = ErrorFullAll(:,:,:,find(ClutterList==-30));
ErrorMaps = reshape(ErrorMaps,3,numel(ll),numel(ll),[],Nalgo);
ErrorMaps = squeeze(mean(ErrorMaps,4));
ErrorMaps = permute(ErrorMaps,[2 3 1 4]);
ErrorMaps = abs(ErrorMaps); %[Nx Nz Npos Nalgo]
% For each algo, display the error map corresponding to the 21x21 grid in z, x, and RMSE.
% The value is averaged with the number of occurence of noise and the number of random
% shift (3x2 per pixel)

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
clb.Position(1) = aa(numel(Algo2Disp)).Position(1) + aa(numel(Algo2Disp)).Position(3)*1.05;clb.Position([2 4]) = aa(numel(Algo2Disp)).Position([2 4]);
clb.Label.String='Absolute error';

for ii = 1:numel(Algo2Disp)
    ialgo = Algo2Disp(ii);
    axes(aa(ii+numel(Algo2Disp)))
    imagesc(ll,ll,ErrorMaps(:,:,1,ialgo)),axis image
    caxis(Caxis.y)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    if ii==1,ylabel('Axial Error');end
end
cbr_z = colorbar;
cbr_z.Position(1) = aa(numel(Algo2Disp)*2).Position(1) + aa(numel(Algo2Disp)*2).Position(3)*1.05;cbr_z.Position([2 4]) = aa(numel(Algo2Disp)*2).Position([2 4]);
cbr_z.Label.String='Absolute error';

for ii = 1:numel(Algo2Disp)
    ialgo = Algo2Disp(ii);
    axes(aa(ii+numel(Algo2Disp)*2))
    imagesc(ll,ll,ErrorMaps(:,:,3,ialgo)),axis image
    caxis(Caxis.RMSE)
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    if ii==1,ylabel('RMSE');end
end
cbr_r = colorbar;
cbr_r.Position(1) = aa(numel(Algo2Disp)*3).Position(1) + aa(numel(Algo2Disp)*3).Position(3)*1.05;cbr_r.Position([2 4]) = aa(numel(Algo2Disp)*3).Position([2 4]);
cbr_r.Label.String='Absolute error';

%% Save fig
figname = 'fig4_ErrorsMaps_30dB';
saveas(fig46,[myfilepath_res filesep figname '.png'])
savefig(fig46,[myfilepath_res filesep figname])

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

m_w = .02;
m_h = .09;
gap_w = .008;
gap_h = .1;

pr_s_w = (1-2*m_w-gap_w*(nx*2-1))/nx/2;
pr_s_hs = .15;
pr_s_h = (1-m_h*2-gap_h-pr_s_hs)/1;

posax = @(ix,ip,isc)[m_w+(ix-1)*pr_s_w+(ix-1)*gap_w,...
    m_h+ip*pr_s_h+isc*pr_s_hs+(ip+isc)*gap_h,...
    pr_s_w,...
    pr_s_hs*(mod(ip+isc,3)==1) + pr_s_h*(1-(mod(ip+isc,3)==1))];

half_fontsize = 9;
f1 = figure(1);clf
f1.Position = [46 600 1479 150];
for ix = 1:nx
    ind =and(listPost2(:,1)==listX(ix),listPost2(:,3)==0);
    ind = find(ind);
    
    axes('Position',posax(ix,1,0))
    plot(listX(ix)*Media.Displacement(1),0,'k.','MarkerSize',10)
    axis image
    xlim([-1 1]*.5);ylim([-1 1]*.5)
    aa=gca;
    aa.XAxisLocation = 'origin';aa.YAxisLocation = 'origin';aa.XAxis.Color = [1 1 1]*.5;aa.YAxis.Color = aa.XAxis.Color;
    if ix==1
        text(-.5,0,'-.5','FontSize',half_fontsize,'VerticalAlignment','middle','HorizontalAlignment','right');
        text(.5,0,'.5','FontSize',half_fontsize,'VerticalAlignment','middle','HorizontalAlignment','left');
        text(0,-.5,'-.5','FontSize',half_fontsize,'VerticalAlignment','top','HorizontalAlignment','center');
        text(0,.5,'.5','FontSize',half_fontsize,'VerticalAlignment','bot','HorizontalAlignment','center');
    end
    % lateral profil
    axes('Position',posax(ix,0,0));hold on
    p0 = plot(rdisp,IQ(rdisp+5,6,ind(1)),'.-','LineWidth',1.2,'MarkerSize',8);
    plot(rdisp,IQ(5,rdisp+6,ind(1)),'LineWidth',p0.LineWidth,'MarkerSize',p0.MarkerSize,'Marker','+')
    grid on;ylim([0 1])
    if ix>0;set(gca,'XTickLabel',{});set(gca,'YTickLabel',{});end
end
%
for ix = 1:nx
    ind =and(listPost2(:,1)==0,listPost2(:,3)==listY(ix));
    ind = find(ind);
    
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
    axes('Position',posax(ix+nx,0,0));hold on
    p0 = plot(rdisp,IQ(rdisp+5,6,ind(1)),'.-','LineWidth',1.2,'MarkerSize',8);
    plot(rdisp,IQ(5,rdisp+6,ind(1)),'LineWidth',p0.LineWidth,'MarkerSize',p0.MarkerSize,'Marker','+')
    grid on;ylim([0 1])
    if ix>0;set(gca,'XTickLabel',{});set(gca,'YTickLabel',{});end
end

%% Save fig
figname = 'fig4_ProfilsIQ3_c';
% saveas(gcf,[myfilepath_res filesep figname '.png'])
% savefig(gcf,[myfilepath_res filesep figname])
print(f1,[myfilepath_res filesep figname],'-dpng','-r300')





