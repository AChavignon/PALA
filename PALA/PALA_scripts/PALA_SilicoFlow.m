%% PALA_SilicoFlow.m : Post Processing - errors and pairings algorithms for IN SILICO FLOW
% Performs ULM on Verasonics simulation. IQ are loaded, then noised at different clutter
% level. For each algorithm, bubbles are detected, localized and tracks. A pairing
% algorithm is used to detect good localization, missed bubbles, and wrong detection.
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
fprintf('Running PALA_SilicoFlow.m\n');t_start=tic;
workingdir = [PALA_data_folder '\PALA_data_InSilicoFlow'];cd(workingdir)
filename = 'PALA_InSilicoFlow';

myfilepath = [workingdir filesep filename];
myfilepath_data = [workingdir filesep 'IQ' filesep filename];

listVar = {'P','PData','Trans','Media','UF','Resource','Receive','filetitle'};
load([myfilepath '_sequence.mat'],'-mat',listVar{:})
ConfigFile = load('PALA_InSilicoFlow_v3_config.mat','MatOut','MyMedia');MatOut_target=ConfigFile.MatOut;

myfilepath_res = [workingdir filesep 'Results'];mkdir(myfilepath_res)
myfilepath_res = [myfilepath_res filesep filename];
%% ULM parameters
ConfigFile.MyMedia.nbPointsMax
NFrames = P.BlocSize*P.numBloc; framerate = P.FrameRate;
res = 10;
ULM = struct('numberOfParticles',40,...% Number of particles per frame. (30-100)
    'res',10,...                    % Resolution factor. Typically 10 for images at lambda/10.
    'max_linking_distance',2,...    % Maximum linking distance between two frames to reject pairing, in pixels units (UF.scale(1)). (2-4 pixel).
    'min_length', 15,...            % Minimum length of the tracks. (5-20)
    'fwhm',[3 3],...                % Size of the mask for localization. (3x3 for pixel at lambda, 5x5 at lambda/2). [fmwhz fmwhx]
    'max_gap_closing', 0,...        % Allowed gap in microbubbles' pairing. (0)
    'size',[PData.Size(1),PData.Size(2),NFrames],...
    'scale',[1 1 1/framerate],...   % Scale [z x t]
    'numberOfFramesProcessed',NFrames,... % Number of processed frames
    'interp_factor',1/res);         % Interpfactor
listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
Nalgo = numel(listAlgo);
%% Simulated Noise parameters
% created a clutter-like noise.
NoiseParam.Power        = -2;   % [dBW]
NoiseParam.Impedance    = .2;   % [ohms]
NoiseParam.SigmaGauss   = 1.5;  % Gaussian filtering
NoiseParam.clutterdB    = -20;  % Clutter level in dB (will be changed later)
NoiseParam.amplCullerdB = 10;   % dB amplitude of clutter

% The code processes the simulated data for various SNR listed above :
ClutterList = [-60 -40 -30 -25 -20 -15 -10]; %dB
%% Display noised IQ example
temp = load([myfilepath_data '_IQ010.mat']);
dB = 20*log10(abs(PALA_AddNoiseInIQ(temp.IQ,NoiseParam))); dB = dB-max(dB(:));

figure(1);imData=imagesc(dB(:,:,1),[-30 0]);colormap gray;axis image;colorbar
for ii=1:100
    imData.CData = dB(:,:,ii);title(ii)
    drawnow;pause(0.01)
end
%% %% Load data and perform detection/localization and tracking algorithms    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--- DETECTION AND LOCALIZATION --- \n\n')
clear Track_tot Track_tot_interp IQ_speckle IQ temp dB
t1=tic;
for iclutter = 1:numel(ClutterList)
    disp(['Computation: SNR at ' num2str(abs(ClutterList(iclutter))) 'dB'])
    NoiseParam.clutterdB = ClutterList(iclutter);
    parfor hhh = 1:P.Nrepeat
        disp(['Bloc ' num2str(hhh) '/' num2str(P.Nrepeat)])
        temp = load([myfilepath_data '_IQ'  num2str(hhh,'%03.0f') '.mat'],'IQ','Media','ListPos');
        Media = temp.Media;ListPos = temp.ListPos;
        IQ = PALA_AddNoiseInIQ(abs(temp.IQ),NoiseParam);

        % save([myfilepath_res '_IQnoise_' num2str(abs(ClutterList(iclutter))) 'dB'],'IQ','ListPos','Media');
        [Track_tot(hhh,:),Track_tot_interp(hhh,:)] = PALA_multiULM(IQ,listAlgo,ULM,PData);
    end
    save([myfilepath_res '_Tracks_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'Track_tot','Track_tot_interp','ULM','P','listAlgo','Nalgo','ClutterList','NoiseParam');disp('Saved')
end
t2=toc(t1);
fprintf('Detection and localization done in %d hours %.1f minutes.\n', floor(t2/60/60), rem(t2/60,60));

%% %% Create MatOut      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--- CREATING MATOUTS --- \n\n')
for iclutter = 1:numel(ClutterList)
    %% load tracks and generate MatOut for each Clutter levels
    disp(['MatOut: SNR at ' num2str(abs(ClutterList(iclutter))) 'dB'])
    MatOut = {};MatOut_vel={};
    load([myfilepath_res '_Tracks_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'Track_tot','Track_tot_interp','ULM','P','listAlgo');disp('Done')

    clear MatOut
    for ialgo = 1:Nalgo
        Track_tot_matout = cat(1,Track_tot_interp{:,ialgo});
        aa = -PData(1).Origin([3 1])+[1 1]*1;  % get origin
        bb = 1./PData(1).PDelta([3 1])*ULM.res;  % fix the size of pixel
        aa(3) = 0;bb(3) = 1;                % for velocity
        Track_matout = cellfun(@(x) (x(:,[1 2 3])+aa).*bb,Track_tot_matout,'UniformOutput',0);
        [MatOut{ialgo},MatOut_vel{ialgo}] = ULM_Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1],'mode','2D_velmean'); %pos in superpix [z x]
        clear Track_tot_matout temp
    end
    MatOut{numel(listAlgo)+1}=MatOut_target;
    % Display MatOut
    figure('Position',[185 150 1550 540],'NumberTitle',1),clf
    a = tight_subplot(2,4,0,0.03,0.0);
    for ii=1:numel(listAlgo)
        axes(a(ii));imagesc((MatOut{ii}).^(1/2));axis image, colormap hot,axis off
        if ii<numel(MatOut),title(listAlgo{ii});else, title('Target');end
        CountMatOut = sum(MatOut{ii}(:));
    end
    linkaxes(a), drawnow
    save([myfilepath_res '_MatOut_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'MatOut','MatOut_vel','ULM','P','PData','listAlgo','Nalgo','CountMatOut','ClutterList');
end

%% %% Pairing Algorithm     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Threshold_pairing = .5; % maximum distance to allow pairings in wavelength
Threshold_TruePos = .25; % maximum distance to accept a True Positive event
fprintf('--- PAIRING ALGORITHM --- \n\n')
for iclutter = 1:numel(ClutterList)
    disp(['Pairing: SNR at ' num2str(abs(ClutterList(iclutter))) 'dB'])
    load([myfilepath_res '_Tracks_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'Track_tot','Track_tot_interp');
    [ErrList,FinalPairs,MissingPoint,WrongLoc,Stat_class] = PALA_Stat_multi(myfilepath_data,Track_tot,listAlgo,PData,Threshold_pairing,Threshold_TruePos);

    save([myfilepath_res '_Stats_multi' num2str(abs(ClutterList(iclutter))) 'dB'],...
        'ErrList','FinalPairs','MissingPoint','WrongLoc','Stat_class','ULM','Threshold_pairing',...
        'P','listAlgo','ClutterList','NoiseParam','Threshold_TruePos');
    disp('Done')
end
t_end = toc(t_start);
fprintf('PALA_SilicoPSF.m performed in %d hours and %.1f minutes\n',floor(t_end/60/60),rem(t_end/60,60));
fprintf('Next, please run PALA_SilicoPSF_fig.m\n')
run('PALA_SilicoFlow_fig.m')
return

%% Additional codes (optional)
[listAlgoName,ListColor,ListMarker,ListShortName] = PALA_GetFormat;
fig43 = figure(20);clf;%fig43.Position = [223 414 1250 340];
postxt_x = -.4;

aa0 = tight_subplot(2,7,[.09 .015],[.07 .01],[.01 .01]);
xedge = linspace(-.6,.6,200);

labelSize = 16;
for ii = 1:numel(listAlgo)
    % for each algo display histograms
    ialgo=ii;axes(aa0(ii));
    axhist_x = histogram((ErrList{ialgo}(:,3)),xedge,'edgecolor','none','Normalization','Probability','FaceColor',[1 1 1]*0.2);
    axhist_x.FaceAlpha = 1;
    axhist_x.FaceColor = ListColor(ii,:);
    ErrMean_x(ialgo) = mean(ErrList{ialgo}(:,3));
    ErrStd_x(ialgo) = std(ErrList{ialgo}(:,3));
    hold on;xlim([-1 1]*xedge(end)),ylim([0 .06])

    text(.02,aa0(ii).YLim(2)-.009,['\sigma=' num2str(round(ErrStd_x(ialgo),2),'%.2f') ],...
        'fontsize',labelSize,'HorizontalAlignment','left','VerticalAlignment','top')
    grid on
    ax = gca;ax.Color = 'none';
    ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';
    ax.YTickLabel = [];ax.XTickLabel = {'-\lambda/4','0','\lambda/4'};ax.XTick = [-.25,0,.25];
    ax.FontSize = labelSize-4;
    ax.TickDir ='in';ax.TickLength = [.04 .2];ax.Box='off';
    ax.GridAlpha = .5; ax.GridLineStyle = '--';
%     title(listAlgoName{ialgo},'color',cfg.title,'fontsize',labelSize,'interpreter','none')
    if ii == 1
        text(postxt_x,mean(ylim),'Lateral error','FontSize',labelSize,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90)
    end
end

for ii= 1:numel(listAlgo)
    ialgo=ii;axes(aa0(ii+7));
    axhist_z = histogram(ErrList{ialgo}(:,2),xedge,'edgecolor','none','Normalization','Probability','FaceColor',[1 1 1]*0.2);
    axhist_z.FaceAlpha = axhist_x.FaceAlpha;
    ErrMean_z(ialgo) = mean(ErrList{ialgo}(:,2));
    ErrStd_z(ialgo) = std(ErrList{ialgo}(:,2));
    axhist_z.FaceColor = ListColor(ii,:);
    hold on
    xlim([-1 1]*xedge(end));ylim([0 .03])
    text(.02,aa0(ii+7).YLim(2)-.001,['\sigma=' num2str(round(ErrStd_z(ialgo),2),'%.2f') ],...
        'fontsize',labelSize,'HorizontalAlignment','left','VerticalAlignment','top')

    grid on
    set(gca, 'yticklabel', []);set(gca, 'xticklabel', {'-\lambda','0','+\lambda'});
    az = gca; az.Color = ax.Color;
    az.XAxisLocation = 'origin';az.YAxisLocation = 'origin';
    az.YTickLabel = [];az.XTickLabel = ax.XTickLabel;az.XTick = ax.XTick;
    az.FontSize = ax.FontSize;az.TickDir = ax.TickDir;az.TickLength = ax.TickLength;az.Box= ax.Box;
    az.GridAlpha = ax.GridAlpha; az.GridLineStyle = ax.GridLineStyle;
    if ii ==1,text(postxt_x,mean(ylim),'Axial error','FontSize',labelSize,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90)
    end
end

%% Display Jaccard Results
idB = 3;
load([myfilepath_res '_Stats_multi' num2str(abs(ClutterList(idB))) 'dB'])
% Stat_class: [Npos_in, Npos_loc, T_pos, F_neg, F_pos, Jaccar] x nb rep

for ialgo=1:numel(listAlgo)
    TP = squeeze(Stat_class(3,ialgo,:));
    FN = squeeze(Stat_class(4,ialgo,:));
    FP = squeeze(Stat_class(5,ialgo,:));
    J_p(ialgo,:) = TP./(FP+TP)*100; % precision
    J_r(ialgo,:) = TP./(FN+TP)*100; % sensitivity
    J_ac(ialgo,:) = TP./(FP+FN+TP)*100; % Jaccard
end

%%
figure(11)
Nalgo =numel(listAlgo);
h1_sensivity = [];p1_sensitivity = [];
for ialgo=1:Nalgo
    for ialgo2=ialgo+1:Nalgo
        [h1_sensivity(ialgo,ialgo2),p1_sensitivity(ialgo,ialgo2),~,stat_ttest] = ttest2(J_r(ialgo,:),J_r(ialgo2,:));
    end
end
h1_jaccard = [];p1_jaccard = [];
for ialgo=1:Nalgo
    for ialgo2=ialgo+1:Nalgo
        [h1_jaccard(ialgo,ialgo2),p1_jaccard(ialgo,ialgo2),~,stat_ttest] = ttest2(J_ac(ialgo,:),J_ac(ialgo2,:));
    end
end
h1_precision = [];p1_precision = [];
for ialgo=1:Nalgo
    for ialgo2=ialgo+1:Nalgo
        [h1_precision(ialgo,ialgo2),p1_precision(ialgo,ialgo2),~,stat_ttest] = ttest2(J_p(ialgo,:),J_p(ialgo2,:));
    end
end
disp('pvalue for sensitivity');disp(p1_sensitivity)
disp('pvalue for Jaccard');disp(p1_jaccard)
disp('pvalue for precision');disp(p1_precision)

%%
figure(56),clf;
subplot 311;hold on
for ialgo=1:Nalgo
    boxplot(J_p(ialgo,:),'positions',ialgo,'labels',listAlgo(ialgo))
end
ca = gca;
ca.XAxis.TickValues = 1:7;ca.XAxis.TickLabels = listAlgo;
xlim([0 8]),ylim([0 100]),title('Precision')

subplot 312;hold on
for ialgo=1:Nalgo
    boxplot(J_r(ialgo,:),'positions',ialgo,'labels',listAlgo(ialgo))
end
ca = gca;
ca.XAxis.TickValues = 1:7;ca.XAxis.TickLabels = listAlgo;
xlim([0 8]),ylim([0 100]),title('Sensitivity')

subplot 313;hold on
for ialgo=1:Nalgo
    boxplot(J_ac(ialgo,:),'positions',ialgo,'labels',listAlgo(ialgo))
end
ca = gca;
ca.XAxis.TickValues = 1:7;ca.XAxis.TickLabels = listAlgo;
xlim([0 8]),ylim([0 100]),title('Jaccard')

