%% ULM PARLA: Post Processing - filtering, localization and tracking for multi algorithms for IN VIVO BRAIN
% Perfoms ULM on rat brain data.
% IQ are loaded, filtered and processed with localization algorithms.
% For each algorithm, bubbles are detected, localised and tracks.
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
SimppleTracker_folder = 'F:\ArthurC\CHAPRO_local\_SIMULATIONS\FightClub_PostPro\SimpleTracker'; % location of the SimpleTracker
PARLA_data_folder = 'F:\ArthurC\Data\Fight_Club\Data'; % path of data
addpath(genpath(PARLA_addons_folder))
addpath(genpath(SimppleTracker_folder))
%% Selected data file and saving folders
workingdir = [PARLA_data_folder '\PARLA_data_InVivoBrain'];
filename = 'PARLA_InVivoRatBrain_';
cd(workingdir)

myfilepath = [workingdir filesep filename];
mydatapath = [workingdir filesep 'IQ' filesep filename]; % add ' num2str(hhh,'%.3d')
trackspath = [workingdir filesep 'Tracks' filesep filename];
savingpath = [workingdir filesep 'fightclub_round2' filesep filename];

% savingpath = [workingdir filesep 'fc_test' filesep filename];
% trackspath = [workingdir filesep 'fc_test' filesep 'Tracks' filesep filename];

% load([myfilepath '_Sequence_param.mat'])
load([mydatapath num2str(1,'%.3d')],'UF','PData');

%% Adapt parameters
load([mydatapath num2str(1,'%.3d')],'IQ');
NFrames = size(IQ,3);
PData.Origin = [0 PData.Size(2)/2*PData.PDelta(2) 0];
P.FrameRate = UF.FrameRateUF;

%% ULM parameters
res = 10;
ULM = struct('numberOfParticles', 60,...  % Number of particle per frame. (30-100)
    'res',10,... % Resolution factor. Typically 10 for images at lambda/10.
    'seuil',[5 UF.NbFrames],... % svd filtergin
    'max_linking_distance',20,...   % Maximum linking distance between two frames to reject pairing, in super-resolved pixels. (20-40).
    'min_length', 15,...% Minimum allowed length of the tracks. (5-20)
    'fwhm',[1 1]*3,... % Size of the mask for localization. (3x3 pour BF à lambda, 5x5 pour lambda/2). [fmwhz fmwhx fmwhz]
    'max_gap_closing', 0,...%Allowed gap in microbbule pairing. (0)
    'size',[PData.Size(1),PData.Size(2),UF.NbFrames],...
    'scale',[1 1 1],... %sacle [z x y t]
    'numberOfFramesProcessed',UF.NbFrames,... %number of processed frames
    'interp_factor',1/res,... % interfactor
    'Traking', 1,... % tracking 1 or 0
    'Velocimetry', 1,... % tracking 1 or 0
    'Mode','PW');
ULM.butter.CuttofFreq = [50 250];          % Cuttof frequency (Hz) for additional filter. Typically [20 300] at 1kHz.
ULM.butter.samplingFreq = 1000;       % Sampling frequency (Hz)
[but_b,but_a] = butter(2,ULM.butter.CuttofFreq/(ULM.butter.samplingFreq/2),'bandpass');

res = ULM.res;

lx = PData.Origin(1) + [0:PData.Size(2)-1].*PData.PDelta(1);
lz = PData.Origin(3) + [0:PData.Size(1)-1].*PData.PDelta(3);

ULM.algorithm = 'wa'; % will be replace later

listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
listAlgo = {'no_shift','wa','radial'};

irecon = '1';
Nalgo = numel(listAlgo);
%% select SVD filtering Noise
hhh=10;
load([mydatapath num2str(hhh,'%.3d')],'IQ');
bulles = SVDfilter(IQ,ULM.seuil);
bulles = filter(but_b,but_a,bulles,[],3);
bulles(~isfinite(bulles))=0;

figure(1)
dB = 20*log10(abs(bulles)); dB = dB-max(dB(:));
imagesc(lx,lz,dB(:,:,20),[-30 0]),colormap gray
colorbar,axis image
clear temp

%% debug
if 0
IQ_filt = SVDfilter(IQ,ULM.seuil);
IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
IQ_filt(~isfinite(IQ_filt))=0;

tic
[Track_test,Track_test_interp,ProcessingTimeTest(:)] = FC_Superloc_multi2(IQ_filt,listAlgo([2 7]),ULM,PData,UF.FrameRateUF);
toc

Track_matout = Track_test_interp{end};
Track_matout = cellfun(@(x) flip((x(:,[2 1]) - PData(1).Origin([1 3])+[1 1]*1)./PData(1).PDelta([1 3])*ULM.res,2),Track_matout,'UniformOutput',0);
MatOut_i = Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1); %pos in superpix [z x]
figure(1)
imagesc(MatOut_i.^(1/3))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load and localize data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Track_tot Track_tot_interp ProcessingTime
tic
parfor hhh = 1:min(UF.Repeat_sequence,12*10)
    disp(hhh)
    tmp = load([mydatapath num2str(hhh,'%.3d')],'IQ');
    IQ_filt = SVDfilter(tmp.IQ,ULM.seuil);tmp = [];
    IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0;
    
    [Track_tot(hhh,:),Track_tot_interp(hhh,:),ProcessingTime(hhh,:)] = FC_Superloc_multi2(IQ_filt,listAlgo,ULM,PData,UF.FrameRateUF);
    % if no parpool
%     [Track_raw,Track_interp,ProTime] = FC_Superloc_multi2(IQ_filt,listAlgo,ULM,PData,UF.FrameRateUF);
%     save([trackspath 'Tracks' num2str(hhh,'%.3d')],'Track_raw','Track_interp','ProTime','ULM','UF','PData','-v6')

end
Tend = toc
disp('Done')
% save([savingpath '_Tracks_multi' irecon],'Track_tot','PData','Track_tot_interp','ProcessingTime','ULM','P','listAlgo','-v7.3','-nocompression')

%save in different files
for hhh = 1:size(Track_tot,1)
    Track_raw = Track_tot(hhh,:);
    Track_interp = Track_tot_interp(hhh,:);
    ProTime = ProcessingTime(hhh,:);
    save([trackspath 'Tracks' num2str(hhh,'%.3d')],'Track_raw','Track_interp','ProTime','ULM','UF','PData','listAlgo','-v6')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create MatOuts     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each algo, create the MatOut density with interpolated tracks for visual analysis,
% and with non interpolated tracks for aliasing index calculation.
MatOutSat = [];
NbrOfLoc = zeros(Nalgo,1);
MatOut = cell(Nalgo,1);MatOut(:)={0};MatOutNoInterp = MatOut;
   
for hhh=1:UF.Repeat_sequence
    %% Generate MatOut density matrix
    disp(hhh)
    load([trackspath 'Tracks' num2str(hhh,'%.3d')],'Track_raw','Track_interp')
    
    for ialgo = 1:Nalgo
        Track_matout = Track_interp{ialgo};
        Track_matout = cellfun(@(x) flip((x(:,[2 1]) - PData(1).Origin([1 3])+[1 1]*1)./PData(1).PDelta([1 3])*ULM.res,2),Track_matout,'UniformOutput',0);
        MatOut_i = Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1); %pos in superpix [z x]
        clear Track_matout
        MatOut{ialgo} = MatOut{ialgo}+MatOut_i;
        MatOutSat(hhh,ialgo) = nnz(MatOut{ialgo}>0); % calculate saturation curve
        
        Track_matout = Track_raw{ialgo};
        Track_matout = cellfun(@(x) flip((x(:,[2 1]) - PData(1).Origin([1 3])+[1 1]*1)./PData(1).PDelta([1 3])*ULM.res,2),Track_matout,'UniformOutput',0);
        MatOut_i = Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1); %pos in superpix [z x]
        MatOutNoInterp{ialgo} = MatOutNoInterp{ialgo}+MatOut_i;
        Track_count = cat(1,Track_matout{:});
        NbrOfLoc(ialgo) = NbrOfLoc(ialgo)+size(Track_count,1);
    end
end
clear Track_raw Track_interp Track_count Track_matout MatOut_i
% uncomment to save results
% save([savingpath '_MatOut_multi' irecon],'MatOut','MatOutSat','NbrOfLoc','ULM','listAlgo','PData','UF','ProcessingTime');
% save([savingpath '_MatOut_multi_nointerp' irecon],'MatOutNoInterp','MatOutSat','NbrOfLoc','ULM','listAlgo','PData','UF','ProcessingTime');

%%
figure(90),clf
a = tight_subplot(2,ceil(numel(listAlgo)/2));
for ii=1:numel(listAlgo)
    axes(a(ii))
    imagesc(MatOut{ii}.^(1/3))
    axis image, colormap hot
    title(listAlgo{ii},'interpreter','none')
    caxis([0 10])
end
linkaxes(a,'xy')

% saveas(gcf,[savingpath '_multiMatOut_disp .png'])


