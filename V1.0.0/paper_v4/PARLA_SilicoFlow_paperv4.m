%% ULM PARLA: Post Processing - errors and pairings algorithms for IN SILICO FLOW
% Perfoms ULM on Verasonics simulation. IQ are loaded, then noised at different clutter
% level. For each algorithm, bubbles are detected, localised and tracks. A pairing
% algorithm is used to detect good localisation, missed bubbles, and wrong detection.
%
% At the end of the code few draft figures will be created.
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
%% Select IQ and Media folder
workingdir = [PARLA_data_folder '\PARLA_data_InSilicoFlow'];
filename = '2D_flow_1203-180447';
cd(workingdir)

myfilepath = [workingdir filesep filename];
myfilepath_data = [workingdir filesep 'IQ' filesep filename];

listVar = {'P','PData','Trans','Media','UF','Resource','Receive','filetitle'};
load([myfilepath '_sequence.mat'],'-mat',listVar{:})
ConfigFile = load('2Dpaper3_config.mat','MatOut','MyMedia');MatOut_target=ConfigFile.MatOut;

suffixe = '3';

myfilepath_res = [workingdir filesep 'config' suffixe filesep filename];
mkdir([workingdir filesep 'config' suffixe])

%% ULM parameters
ConfigFile.MyMedia.nbPointsMax
NFrames = P.BlocSize*P.numBloc;
res = 10;
ULM = struct('numberOfParticles', 40,...  % Number of particle per frame. (30-100)
    'res',10,... % Resolution factor. Typically 10 for images at lambda/10.
    'seuil',[1 NFrames],... % svd filtergin
    'max_linking_distance',20,...   % Maximum linking distance between two frames to reject pairing, in super-resolved pixels. (20-40).
    'min_length', 15,...% Minimum allowed length of the tracks. (5-20)
    'fwhm',[3 3],... % Size of the mask for localization. (3x3 pour BF à lambda, 5x5 pour lambda/2). [fmwhz fmwhx fmwhz]
    'max_gap_closing', 0,...%Allowed gap in microbbule pairing. (0)
    'size',[PData.Size(1),PData.Size(2),NFrames],...
    'scale',[1 1 1],... %sacle [z x y t]
    'numberOfFramesProcessed',NFrames,... %number of processed frames
    'interp_factor',1/res,... % interfactor
    'Traking', 1,... % tracking 1 or 0
    'Velocimetry', 1,... % tracking 1 or 0
    'Mode','PW');
res = ULM.res;

lx = PData.Origin(1) + [0:PData.Size(2)-1]*PData.PDelta(1);
lz = PData.Origin(3) + [0:PData.Size(1)-1]*PData.PDelta(3);

ULM.algorithm = 'wa'; % SAF radial interp

listAlgo = {'no_localization','SAF','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
listAlgo = {'no_shift','wa','radial'};

% listAlgo = {'no_localization','wa','interp_cubic','radial'};
% listAlgo = {'radial'};

%% Simulated Noise paramters
% created a clutter-like noise.
NoiseParam.Power        = -2;   % [dBW]
NoiseParam.Impedance    = .2;   % [ohms]
NoiseParam.SigmaGauss   = 1.5;  % gaussian filtering
NoiseParam.clutterdB    = -20;  % clutter level in dB (will be changed later)
NoiseParam.amplCullerdB = 10;   % dB amplitude of clutter

% The code processes the simulated data for various SNR listed above :
ClutterList = [-60 -40 -30 -25 -20 -15 -10]; %dB

%% Display noised IQ example
temp = load([myfilepath_data '_IQ'  num2str(10,'%03.0f') '.mat']);

IQ = abs(temp.IQ);
% IQ_speckle = IQ+imgaussfilt(max(IQ(:))*10^(NoiseParam.clutterdB/20)+reshape(wgn(numel(IQ),1,NoiseParam.Power,NoiseParam.Impedance),size(IQ,1),size(IQ,2),[])*max(IQ(:))*10^((NoiseParam.amplCullerdB+NoiseParam.clutterdB)/20),NoiseParam.SigmaGauss);
IQ_speckle = AddNoiseInIQ(IQ,NoiseParam);

figure(1)
dB = 20*log10(abs(IQ_speckle(:,:,10))); dB = dB-max(dB(:));
imagesc(lx,lz,dB,[-30 0]),colormap gray
colorbar,axis image

%% Display noise
figure(1);
dB = 20*log10(abs(IQ_speckle)); dB = dB-max(dB(:));
imData=imagesc(lx,lz,dB(:,:,1),[-30 0]);colormap gray;axis image
colorbar

for ii=1:size(IQ_speckle,3)
    imData.CData = dB(:,:,ii);
    title(ii)
    drawnow;pause(0.01)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data dans performs detection/localization and tracking algorithms    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Track_tot Track_tot_interp IQ_speckle IQ temp
hhh=1;
parpool(20)
for iclutter = 1:numel(ClutterList)
    NoiseParam.clutterdB = ClutterList(iclutter);
    parfor hhh = 1:min(P.Nrepeat,999)
        disp(hhh)
        temp = load([myfilepath_data '_IQ'  num2str(hhh,'%03.0f') '.mat'],'IQ','Media','ListPos');
        Media = temp.Media;ListPos = temp.ListPos;
        IQ = AddNoiseInIQ(abs(temp.IQ),NoiseParam);

        % save([myfilepath_res '_IQnoise_' num2str(abs(ClutterList(iclutter))) 'dB'],'IQ','ListPos','Media');
        [Track_tot(hhh,:),Track_tot_interp(hhh,:)] = FC_Superloc_multi2(IQ,listAlgo,ULM,PData,P.FrameRate);
    end
    disp('Done')
    save([myfilepath_res '_Tracks_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'Track_tot','Track_tot_interp','ULM','P','listAlgo','ClutterList','NoiseParam');disp('Done')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create MatOut      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iclutter = 1:numel(ClutterList)
    %% load tracks and generate MatOut for each Clutter levels
    disp(iclutter)
    load([myfilepath_res '_Tracks_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'Track_tot','Track_tot_interp','ULM','P','listAlgo');disp('Done')

    clear MatOut
    for iaglo = 1:numel(listAlgo)
        Track_tot_matout = cat(1,Track_tot_interp{:,iaglo});
        Track_tot_matout = cellfun(@(x) flip((x(:,[2 1]) - PData(1).Origin([1 3])+[1 1]*1)./PData(1).PDelta([1 3])*ULM.res,2),Track_tot_matout,'UniformOutput',0);
        MatOut{iaglo} = Track2MatOut(Track_tot_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1); %pos in superpix [z x]
        clear Track_tot_matout temp
    end
    
    %% Display MatOut
    MatOut{numel(listAlgo)+1}=MatOut_target;
    figure('Position',[185 72 1320 906]),clf
    a = tight_subplot(3,3,0,0.03,0.0);
    for ii=1:numel(listAlgo)
        axes(a(ii))
        imagesc((MatOut{ii}).^(1/2))
        axis image, colormap hot,axis off
        if ii<numel(MatOut),title(listAlgo{ii});else title('Target');end
        CountMatOut = sum(MatOut{ii}(:));
    end
    linkaxes(a)
    saveas(gcf,[myfilepath_res '_multiMatOut' num2str(abs(ClutterList(iclutter))) 'dB.png'])
    
    save([myfilepath_res '_MatOut_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'MatOut','ULM','P','listAlgo','CountMatOut','ClutterList');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pairing Algorithm     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Threshold_pairing = .25; % maximum distance to allow pairings

for iclutter = 1:numel(ClutterList)
    clear ErrList Npos_in Npos_loc T_pos F_neg F_pos
    clear ErrList FinalPairs MissingPoint WrongLoc Stat_class
    load([myfilepath_res '_Tracks_multi_' num2str(abs(ClutterList(iclutter))) 'dB'],'Track_tot','Track_tot_interp');disp('Done')
    [ErrList,FinalPairs,MissingPoint,WrongLoc,Stat_class] = FC_Stat_multi(myfilepath_data,Track_tot,listAlgo,PData,Threshold_pairing);
    save([myfilepath_res '_Stats_multi' num2str(abs(ClutterList(iclutter))) 'dB'],'ErrList','FinalPairs','MissingPoint','WrongLoc','Stat_class','ULM','Threshold_pairing','P','listAlgo','ClutterList','NoiseParam');disp('Done')
end

%% Display Loc Errors
% ErrList
idB = 1;
load([myfilepath_res '_Stats_multi' num2str(abs(ClutterList(idB))) 'dB'])

for ialgo=1:numel(listAlgo)
    figure(46+ialgo);clf
        
    subplot(3,2,1)
    h = histogram(ErrList{ialgo}(:,1),200,'edgecolor','none','Normalization','Probability');
    xlim([0 1])
    
    title(['RMSE dist ' num2str(mean(ErrList{ialgo}(:,1))) '\lambda'])
    grid on,xlabel('RMSE [\lambda]')
    
    subplot(3,2,2)
    cdfplot(ErrList{ialgo}(:,1))
    hold on
    
    title('RMSE cdf')
    grid on,xlabel('RMSE [\lambda]')
    
    subplot(3,2,3)
    histogram(ErrList{ialgo}(:,2),200,'edgecolor','none','Normalization','Probability');
    title('Err z')
    grid on,xlabel('Err_z [\lambda]')
    xlim([-1 1])
    
    subplot(3,2,4)
    histogram(ErrList{ialgo}(:,3),200,'edgecolor','none','Normalization','Probability');
    title('Err x')
    grid on,xlabel('Err_x [\lambda]')
    xlim([-1 1])
    
    subplot 313,hold on
    myBoxPlotFC({'RMSE','z','x'},ErrList{ialgo})
    ylabel('\lambda','Rotation',0),grid on,grid minor
    title('Localisation errors')
    ylim([-.3 .6])
    
    suptitle([listAlgo{ialgo} ' - Jaccard indice ' num2str(round(100*Stat_class(6,ialgo),1))])
%     saveas(gcf,[myfilepath_res '_ErrorHist_' listAlgo{ialgo} '_' num2str(abs(ClutterList(idB))) 'dB.png'])
end

%% Compare error and RMSE
figure(12),clf
% error z
subplot 312
dd = cellfun(@(x) (x(:,3)-mean(x(:,3),1)),ErrList,'UniformOutPut',false);
myBoxPlotFC(listAlgo,dd)
ylabel('\lambda','Rotation',0),grid on,grid minor,drawnow
title('Error Z')

% error x
subplot 311
dd = cellfun(@(x) (x(:,2)-mean(x(:,2),1)),ErrList,'UniformOutPut',false);

myBoxPlotFC(listAlgo,dd)
ylabel('\lambda','Rotation',0),grid on,grid minor,drawnow
title('Error X')

% error RMSE
subplot 313
dd = cellfun(@(x) x(:,1),ErrList,'UniformOutPut',false);
myBoxPlotFC(listAlgo,dd)
%      myBoxPlotFC_multi(listAlgo,dd,10,hot',12,1)
ylabel('\lambda','Rotation',0),grid on,grid minor,drawnow
title('Root Mean Square Error')

% saveas(gcf,[myfilepath_res '_ErrorRMSE_' listAlgo{ialgo} '_' num2str(abs(ClutterList(idB))) 'dB.png'])

%% Display Jaccard Results
figure(75)
% [Npos_in x Npos_loc x T_pos x F_neg x F_pos x Jaccar]
subplot 221
bar(Stat_class(3,:))
title('T Pos')

subplot 222
bar(Stat_class(4,:))
title('F neg')

subplot 223
bar(Stat_class(5,:))
title('F pos')

subplot 224
bar(Stat_class(2,:))
title('Npos loc')

for ialgo=1:numel(listAlgo)
TP = Stat_class(3,ialgo);
FN = Stat_class(4,ialgo);
FP = Stat_class(5,ialgo);

J_p(ialgo) = TP/(FP+TP)*100; % precision
J_r(ialgo) = TP/(FN+TP)*100; % sensitivity
J_ac(ialgo) = TP/(FP+FN+TP)*100; % Jaccard
end

%% Display for 1 frame
hhh = 1;
ifr = 10;

temp1 = load([myfilepath_data '_IQ'  num2str(hhh,'%03.0f') '.mat'],'ListPos','IQ');

figure(2);clf
db = 20*log10(abs(temp1.IQ(:,:,ifr)));
imagesc(lx,lz,db-max(db(:)),[-40 0])
% cm = flip(gray);
cm = gray;
colormap(cm);colorbar
hold on

ialgo = 3;
msize = 5;

plot(temp1.ListPos(FinalPairs{ialgo}{ifr}(:,1),1,ifr),temp1.ListPos(FinalPairs{ialgo}{ifr}(:,1),3,ifr),...
    'b.','markersize',msize,'DisplayName','True Positive');

TT = cell2mat(Track_tot{hhh,ialgo});
TT = TT(TT(:,3)==ifr,:);

plot(TT(FinalPairs{ialgo}{ifr}(:,2),2),TT(FinalPairs{ialgo}{ifr}(:,2),1),...
    'bx','markersize',msize,'DisplayName','True Positive loc');

plot(temp1.ListPos(MissingPoint{ialgo}{ifr},1,ifr),temp1.ListPos(MissingPoint{ialgo}{ifr},3,ifr),...
    'r.','markersize',msize,'DisplayName','True negative (missing points)')

plot(TT(WrongLoc{ialgo}{ifr},2),TT(WrongLoc{ialgo}{ifr},1),...
    'rd','markersize',msize,'DisplayName','False negatieve')

title([num2str(ifr) ' algo ' listAlgo{ialgo}],'interpreter','none')
xlabel('\lambda');ylabel('\lambda')
grid on
legend('Location','southoutside')

axis image
xlim([-68 68])

% saveas(gcf,[workingdir filesep 'img_4' filesep filename '_DisplayLoc' suffixe '.png'])

