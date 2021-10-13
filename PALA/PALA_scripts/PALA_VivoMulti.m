%% PALA_VivoMulti.m : Post Processing - filtering, localization and tracking for multi algorithms for IN VIVO DATASETS
% Performs ULM on the selected dataset. IQ are loaded, filtered and processed
% with localization algorithms. For each algorithm, bubbles are detected, localised and tracks.
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
[listAlgoName,ListColor,ListMarker,ListShortName] = PALA_GetFormat;
% SELECT DATASET NUMBER
DataSetNumber = 1; % 1 2 or 3

%%
workingdir = [PALA_data_folder '\'];
% ULMparam is a summmary of filtering and processing parameters
% ULMparam = [nb_part,SVD_cutoff,max_linking_distance,min_length,butter_1,butter_end,NLocalMax]
switch DataSetNumber
    case 1 % Rat brain with bolus
        workingdir = [workingdir '\PALA_data_InVivoRatBrainBolus'];
        filename = 'PALA_InVivoRatBrainBolus_';
        ULMparam = [100,10,2,10,50,250,4];
    case 2 % Mouse Tumor
        workingdir = [workingdir '\PALA_data_InVivoMouseTumor'];
        filename = 'PALA_InVivoMouseTumor_';
        ULMparam = [200,10,1,8,25,249,4];
    case 3 % Rat kidney
        workingdir = [workingdir '\PALA_data_InVivoRatKidney'];
        filename = 'PALA_InVivoRatKidney_';
        ULMparam = [100,10,1.5,15,30,250,3];
    otherwise
        error('Select DataSetNumber = {1,2,3}.')
end
cd(workingdir)

%% Selected data file and saving folders
fprintf('--- Running PALA_VivoMulti.m with dataset:  %s ---\n\n',filename);t_start =tic;
myfilepath = [workingdir filesep filename];
mydatapath = [workingdir filesep 'IQ' filesep filename]; % add ' num2str(hhh,'%.3d')
trackspath = [workingdir filesep 'Tracks' filesep filename];mkdir(fileparts(trackspath))
savingpath = [workingdir filesep 'Results' filesep filename];mkdir(fileparts(savingpath))

IQfiles = dir([mydatapath '*.mat']);Nbuffers = numel(IQfiles);
load([IQfiles(1).folder filesep IQfiles(1).name],'UF');

%% Adapt parameters
load([IQfiles(1).folder filesep IQfiles(1).name],'IQ','PData');
NFrames = size(IQ,3);
PData.Size = size(IQ);
PData.PDelta = [1 1 1]; % isotropic pixels
PData.Origin = [0 PData.Size(2)/2*PData.PDelta(2) 0];
framerate = UF.FrameRateUF;

%% ULM parameters
res = 10;
ULM = struct('numberOfParticles', ULMparam(1),...% Number of particles per frame. (30-100)
    'res',10,...                        % Resolution factor. Typically 10 for images at lambda/10.
    'SVD_cutoff',[ULMparam(2) NFrames],...   % svd filtering
    'max_linking_distance',ULMparam(3),... % Maximum linking distance between two frames to reject pairing, in pixels units (UF.scale(1)). (2-4 pixel).
    'min_length', ULMparam(4),...       % Minimum length of the tracks. (5-20)
    'fwhm',[1 1]*3,...                  % Size of the mask for localization. (3x3 for pixel at lambda, 5x5 at lambda/2). [fmwhz fmwhx]
    'max_gap_closing', 0,...            % Allowed gap in microbubbles' pairing. (0)
    'size',[PData.Size(1),PData.Size(2),NFrames],...
    'scale',[1 1 1/framerate],...     % Scale [z x t]
    'numberOfFramesProcessed',NFrames,...% number of processed frames
    'interp_factor',1/res);             % interpfactor

ULM.butter.CuttofFreq = [ULMparam(5) ULMparam(6)];% Cuttof frequency (Hz) for additional filter. Typically [20 300] at 1kHz.
ULM.butter.samplingFreq = framerate;         % Sampling frequency (Hz)
[but_b,but_a] = butter(2,ULM.butter.CuttofFreq/(ULM.butter.samplingFreq/2),'bandpass');
ULM.parameters.NLocalMax = ULMparam(7);

lx = PData.Origin(1) + [0:PData.Size(2)-1].*PData.PDelta(1);
lz = PData.Origin(3) + [0:PData.Size(1)-1].*PData.PDelta(3);

listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
Nalgo = numel(listAlgo);
%% select SVD filtering Noise
hhh=min(Nbuffers,19);
load([IQfiles(hhh).folder filesep IQfiles(hhh).name],'IQ');
IQ_filt = SVDfilter(IQ,ULM.SVD_cutoff);
IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
IQ_filt(~isfinite(IQ_filt))=0;

BullesSNR = abs(IQ_filt(:,:,10));
LocalMax = imregionalmax(BullesSNR);ValMax = sort(BullesSNR(LocalMax),'descend');
SNRmean = 20*log10(mean(ValMax(1:10))/mean(BullesSNR,'all'));clear BullesSNR LocalMax ValMax

figure(1)
dB = 20*log10(abs(IQ_filt)); dB = dB-max(dB(:));
imagesc(lx,lz,dB(:,:,20),[-30 0]),colormap gray
colorbar,axis image;clear temp

%% %% Load and localize data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Track_tot Track_tot_interp ProcessingTime IQ_filt IQ dB
fprintf('--- ULM PROCESSING --- \n\n');t1=tic;
parfor hhh = 1:min(Nbuffers,999)
    fprintf('Processing bloc %d/%d\n',hhh,Nbuffers);
    tmp = load([IQfiles(hhh).folder filesep IQfiles(hhh).name],'IQ');
    tmp.IQ(:,[1 end])=0;tmp.IQ([1 end],:)=0;
    IQ_filt = SVDfilter(tmp.IQ,ULM.SVD_cutoff);tmp = [];
    IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0;
    
    %     [Track_tot(hhh,:),Track_tot_interp(hhh,:),ProcessingTime(hhh,:)] = PALA_multiULM(IQ_filt,listAlgo,ULM,PData);
    % Data will be written in a '.mat' file to free RAM
    [~,~] = PALA_multiULM(IQ_filt,listAlgo,ULM,PData,'savingfilename',[trackspath 'Tracks' num2str(hhh,'%.3d') '.mat']);
end
Tend = toc(t1);
fprintf('ULM done in %d hours %.1f minutes. (for all localization algorithms) \n', floor(Tend/60/60), rem(Tend/60,60));
save([trackspath 'Tracks' num2str(1,'%.3d') '.mat'],'Tend','SNRmean','-append');
% save([savingpath '_Tracks_multi'],'Track_tot','PData','Track_tot_interp','ProcessingTime','ULM','P','listAlgo','-v7.3','-nocompression')
clear tmp IQ_filt

%% %% Create MatOuts     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--- CREATING MATOUTS --- \n\n')
% for each algorithm, creates the MatOut density with interpolated tracks for visual analysis,
% and with non interpolated tracks for aliasing index calculation.
MatOutSat = [];
NbrOfLoc = zeros(Nalgo,1);ProcessingTime = zeros(Nbuffers,Nalgo);
MatOut = cell(Nalgo,1);MatOut(:)={0};MatOutNoInterp = MatOut;MatOut_vel = MatOut;

hwait = waitbar(0,'Bluiding intensity renderings','Name','Building matouts');
for hhh=1:min(Nbuffers,999) % Generate MatOut density matrix
    waitbar(hhh/Nbuffers,hwait);
    load([trackspath 'Tracks' num2str(hhh,'%.3d')],'Track_raw','Track_interp','ProTime')
    aa = -PData(1).Origin([3 1])+[1 1]*1;  %get the good origin
    bb = [1./PData(1).PDelta([3 1])*ULM.res];  %fix the size of pixel
    aa(3) = 0;bb(3) = 1; % for velocity
    for ialgo = 1:Nalgo
        % MatOut and MatOutVel rendering with interpolated tracks
        Track_matout = Track_interp{ialgo};
        Track_matout = cellfun(@(x) (x(:,[1 2 3])+aa).*bb,Track_matout,'UniformOutput',0);
        [MatOut_i,MatOut_vel_i] = ULM_Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1,'mode','2D_velmean'); %pos in superpix [z x]
        clear Track_matout
        MatOut_vel{ialgo} = MatOut_vel{ialgo}.*MatOut{ialgo}+MatOut_vel_i.*MatOut_i; % weighted summation
        MatOut{ialgo} = MatOut{ialgo}+MatOut_i;
        MatOut_vel{ialgo}(MatOut{ialgo}>0) = MatOut_vel{ialgo}(MatOut{ialgo}>0)./MatOut{ialgo}(MatOut{ialgo}>0); % average velocity
        MatOutSat(hhh,ialgo) = nnz(MatOut{ialgo}>0); % compute saturation curve
        
        % MatOut without interpolation, for gridding index
        Track_matout = Track_raw{ialgo};
        Track_matout = cellfun(@(x) (x(:,[1 2])+aa(1:2)).*bb(1:2),Track_matout,'UniformOutput',0);
        MatOut_i = ULM_Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1); %pos in superpix [z x]
        MatOutNoInterp{ialgo} = MatOutNoInterp{ialgo}+MatOut_i;
        Track_count = cat(1,Track_matout{:});clear Track_matout
        NbrOfLoc(ialgo) = NbrOfLoc(ialgo)+size(Track_count,1);
        
    end
    ProcessingTime(hhh,:) = ProTime;
end
clear Track_raw Track_interp Track_count Track_matout MatOut_i
save([savingpath 'MatOut_multi'],'MatOut','MatOut_vel','MatOutSat','NbrOfLoc','ULM','listAlgo','PData','Nalgo','UF');
save([savingpath 'MatOut_multi_nointerp'],'MatOutNoInterp','MatOutSat','NbrOfLoc','ULM','listAlgo','PData','UF','Nalgo','ProcessingTime');
close(hwait);clear hwait

%% %% Displaying %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([savingpath 'MatOut_multi']);load([savingpath 'MatOut_multi_nointerp']);

switch DataSetNumber
    case 1 % Rat brain with bolus
        MatZoom = [582 775 432 662]; %zoom start_x end_x start_y end_y
        IntPower = 1/3;SigmaGauss=0;
    case 2 % Mouse Tumor
        MatZoom = [373 590 735 930];
        IntPower = 1/2;SigmaGauss=.7;
    case 3 % Rat kidney
        MatZoom = [510 760 630 870];
        IntPower = 1/2;SigmaGauss=.7;
    otherwise
        error('Select DataSetNumber = {1,2,3}.')
end

figure(90),clf
a = tight_subplot(2,ceil(Nalgo/2),.01,[.005 .02],.005);
for ialgo=1:Nalgo
    axes(a(ialgo))
    im=imagesc(MatOut{ialgo}.^IntPower);
    if SigmaGauss>0,im.CData = imgaussfilt(im.CData,SigmaGauss);end
    axis image, colormap hot,axis off
    title(listAlgoName{ialgo},'interpreter','none')
    caxis(caxis*.8)
end
suptitle(filename)
linkprop(a,{'CLim','XLim','YLim'})
print([savingpath 'multiMatOut_disp'],'-dpng','-r400')

%% Display all MatOut intensity
wv = 1540/UF.TxFreq*1e-3; %[mm]
figure(91),clf
a = tight_subplot(2,ceil(Nalgo/2));
for ialgo=1:Nalgo
    axes(a(ialgo))
    imagesc(MatOut_vel{ialgo}.*wv) % in [mm/s]
    axis image, colormap jet,axis off
    title(listAlgoName{ialgo},'interpreter','none')
    caxis([0 40])
end
suptitle(filename)
linkprop(a,{'CLim','XLim','YLim'})
print([savingpath 'multiMatVel_disp'],'-dpng','-r400')

%% Scores
% Saturation Value
clear Score
Score.Saturation  = MatOutSat(end,:)/numel(MatOut{1});
% Processing Time
ProTime = sum(ProcessingTime,1);
Score.Time = 1-log10(ProTime./min(ProTime))/2;

% Gridding index
Gridding_index = PALA_ComputeGriddingIndex(MatOutNoInterp,1);
Score.Gridding = 1-Gridding_index/30; % [100 at 0; 0 at -30]

%% %% Displaying - Radial Full Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelSize = 12;
wv = 1540/UF.TxFreq*1e-3; % wavelength in mm;
ScaleBar = 100*PData.PDelta(1)/ULM.res*wv;

fig38 = figure(38);clf;
ialgo = 7; % select radial algorithm
im=imagesc(imgaussfilt(MatOut{ialgo}.^IntPower,.01));
if SigmaGauss>0,im.CData = imgaussfilt(im.CData,SigmaGauss);end
aa = gca;axis image,axis off, colormap(hot(128))
posScale(1) = 60; posScale(2) = size(MatOut{ialgo},1)-70;
tscale = text(posScale(1)+50, posScale(2),'1 mm','color','w','fontsize',labelSize-2,'VerticalAlignment','bot','HorizontalAlignment','center');
aa.CLim=aa.CLim.*0.8;MatClim = aa.CLim;

im.CData(posScale(2)+[-2:2],posScale(1)+[1:100])=aa.CLim(2);
im.CData = min(im.CData,aa.CLim(2));
im.CData(MatZoom(1):MatZoom(2),MatZoom(3)+[0:2])=aa.CLim(2);im.CData(MatZoom(1):MatZoom(2),MatZoom(4)-[0:2])=aa.CLim(2);
im.CData(MatZoom(1)+[0:2],MatZoom(3):MatZoom(4))=aa.CLim(2);im.CData(MatZoom(2)-[0:2],MatZoom(3):MatZoom(4))=aa.CLim(2);
title([filename ': localization algorithm ' listAlgo{ialgo}],'Interpreter','none')

figname = 'fig7_matout_radial_full';
saveas(fig38,[savingpath figname '.png'])
WriteTif(im.CData,aa.Colormap,[savingpath figname '.tif'],'caxis',MatClim,'Overwrite',1)

%% Save all
for ialgo=1:Nalgo
    im=imagesc(imgaussfilt(MatOut{ialgo}.^IntPower,.01));
    if SigmaGauss>0,im.CData = imgaussfilt(im.CData,SigmaGauss);end
    aa = gca;axis image,axis off, colormap(hot(128))
    posScale(1) = 60; posScale(2) = size(MatOut{ialgo},1)-70;
    aa.CLim=aa.CLim.*0.8;
    im.CData(posScale(2)+[-2:4],posScale(1)+[1:100])=aa.CLim(2);
    im.CData = min(im.CData,aa.CLim(2));
    figname = ['figsup_' filename ListShortName{ialgo}];
    
    MatText = zeros(60,110);
    MatText = insertText(MatText,[-20 -30],ListShortName{ialgo},'TextColor','white','BoxOpacity',0,'FontSize',70,'AnchorPoint','LeftTop');
    MatText = (sum(MatText,3)>.5).*aa.CLim(2);
    im.CData([1:size(MatText,1)],[1:size(MatText,2)])=MatText;
    figure(38);WriteTif(im.CData,aa.Colormap,[savingpath figname '.tif'],'caxis',caxis,'Overwrite',1)
end

%%
fig39 = figure(39);clf;fig39.Units = 'centimeters';fig39.Position = [7 2 8 1];
caxis(MatClim);colormap(aa.Colormap)
clb = colorbar('Location','south','Position',[.04 .01 .9 .3]);
val = [0:10:100 200:100:1000];
clb.Ticks = val.^IntPower;
clb.TickLabels = val;clb.TickLabels([3:5 7:10 13:14 16:19],:)=' ';
clb.AxisLocation = 'in';clb.Label.String='Counts';clb.TickLength = .02;
clb.Label.Position = [.9 1.5];clb.FontSize = labelSize-4;axis off
figname = 'fig7_matout_radial_full_clb';drawnow
print([savingpath figname],'-dpng','-r400')
print(fig39,[savingpath figname '.eps'],'-depsc2','-tiff','-r300','-painters')

%% Zoom
fig40 = figure(40);clf;fig40.Position = [5 696 1642 300];
a = tight_subplot(1,5,[0.01 0.005],0.00,0.00);
listAlgoMatOut = [2 3 5 6 7];
for ii=1:numel(listAlgoMatOut)
    ialgo=listAlgoMatOut(ii);axes(a(ii))
    im=imagesc(MatOut{ialgo}(MatZoom(1):MatZoom(2),MatZoom(3):MatZoom(4)).^IntPower);
    if SigmaGauss>0,im.CData = imgaussfilt(im.CData,SigmaGauss);end
    axis image, colormap(hot(128));a(ii).CLim = MatClim;axis off
    im.CData = min(im.CData,MatClim(2));
    
    if ii==1
        posScale(1) = 10; posScale(2) = round(a(ii).YLim(2)-10);hold on
        im.CData(posScale(2)+[0:2],posScale(1)+[1:100])=MatClim(2);
    end
    title(listAlgoName{ialgo})
    figname = ['fig7_matout_zoom_' ListShortName{ialgo}];drawnow
    WriteTif(im.CData,hot,[savingpath figname '.tif'],'caxis',caxis,'overwrite',1)
end
linkaxes(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display score    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ListScore = [Score.Time(:),Score.Gridding(:),Score.Saturation(:)];
lgdname = {'Time','Gridding','Saturation'};
lgdshort = {'Time','Grid','Sat'};

Weight_coef = [1 1 1];
ListScore = ListScore.*Weight_coef./sum(Weight_coef)*100;
GlobalScore = sum(ListScore,2);

fig02 = figure(2);clf;fig02.Position(3:4) = [360 210];axes('Position',[.01 .1 .98 .82])

bb = barh(ListScore,1,'stacked');
for ii=1:numel(bb)
    bb(ii).FaceColor='flat';
    bb(ii).CData = ListColor(:,:)*1;
    bb(ii).FaceAlpha = 1-.7*((ii-1)/(numel(bb)-1));
end

ylim([.5 9]);xlim([0 100]);grid on
aa = gca;set(gca,'ytick',[]);aa.XAxis.TickDirection='both';aa.Color ='none';
aa.Box = 'off';aa.XAxis.TickValues=[0:20:100];

% text(GlobalScore+3,[1:7],ListShortName,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize)
ScoreString = arrayfun(@(x) num2str(round(GlobalScore(x),1)),[1:Nalgo],'UniformOutput',false);
text(GlobalScore+1,[1:7],ScoreString,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',labelSize)

lgd =legend(bb,lgdname,'Location','northwest','Orientation','horizontal');
txtpos = [0 cumsum(ListScore(end,1:end-1))] + ListScore(end,:)/2;
text(txtpos,7*ones(1,numel(txtpos)),lgdshort,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8)

for ii=1:numel(lgd.String),lgd.String{ii} = [lgd.String{ii}];end
title(filename,'Interpreter','none')
figname = ['fig7_score'];drawnow
print([savingpath figname ],'-dpng','-r300')
print(fig02,[savingpath figname '.eps'],'-depsc2','-tiff','-r300','-painters')

t_end = toc(t_start);
fprintf('PALA_VivoMulti.m performed in %d hours and %.1f minutes\n',floor(t_end/60/60),rem(t_end/60,60));
