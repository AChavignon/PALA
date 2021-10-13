%% PALA_VivoBrain.m : Post Processing - filtering, localization and tracking for multi algorithms for IN VIVO BRAIN
% Performs ULM on rat brain data.
% IQ are loaded, filtered and processed with localization algorithms.
% For each algorithm, bubbles are detected, localized and tracks.
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
%% Selected data file and saving folders
fprintf('Running PALA_VivoBrain.m\n');t_start =tic;
workingdir = [PALA_data_folder '\PALA_data_InVivoRatBrain'];cd(workingdir)
filename = 'PALA_InVivoRatBrain_';

myfilepath = [workingdir filesep filename];
mydatapath = [workingdir filesep 'IQ' filesep filename]; % add ' num2str(hhh,'%.3d')
trackspath = [workingdir filesep 'Tracks' filesep filename]; mkdir(fileparts(trackspath))
savingpath = [workingdir filesep 'Results' filesep filename]; mkdir(fileparts(savingpath))

% load([myfilepath '_Sequence_param.mat'])
IQfiles = dir([mydatapath '*.mat']);
Nbuffers = numel(IQfiles);
load([IQfiles(1).folder filesep IQfiles(1).name],'UF','PData');

%% Adapt parameters
load([IQfiles(1).folder filesep IQfiles(1).name],'IQ');
NFrames = size(IQ,3);
PData.Origin = [0 PData.Size(2)/2*PData.PDelta(2) 0];
framerate = UF.FrameRateUF;

%% ULM parameters
res = 10;
ULM = struct('numberOfParticles', 70,...% Number of particles per frame. (30-100)
    'res',10,...                        % Resolution factor. Typically 10 for images at lambda/10.
    'SVD_cutoff',[5 UF.NbFrames],...    % svd filtering
    'max_linking_distance',2,...        % Maximum linking distance between two frames to reject pairing, in pixels units (UF.scale(1)). (2-4 pixel).
    'min_length', 15,...                % Minimum length of the tracks. (5-20)
    'fwhm',[1 1]*3,...                  % Size of the mask for localization. (3x3 for pixel at lambda, 5x5 at lambda/2). [fmwhz fmwhx]
    'max_gap_closing', 0,...            % Allowed gap in microbubbles pairing. (0)
    'size',[PData.Size(1),PData.Size(2),UF.NbFrames],...
    'scale',[1 1 1/framerate],...       % Scale [z x t]
    'numberOfFramesProcessed',UF.NbFrames,... % Number of processed frames
    'interp_factor',1/res...            % interpfactor
    );
ULM.butter.CuttofFreq = [50 250];       % Cut off frequency (Hz) for additional filter. Typically [20 300] at 1kHz.
ULM.butter.samplingFreq = framerate;         % Sampling frequency (Hz)
[but_b,but_a] = butter(2,ULM.butter.CuttofFreq/(ULM.butter.samplingFreq/2),'bandpass');
ULM.parameters.NLocalMax = 3;           % Safeguard on the number of maxLocal in the fwhm*fwhm grid (3 for fwhm=3, 7 for fwhm=5)
res = ULM.res;

lx = PData.Origin(1) + [0:PData.Size(2)-1].*PData.PDelta(1);
lz = PData.Origin(3) + [0:PData.Size(1)-1].*PData.PDelta(3);
listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
Nalgo = numel(listAlgo);
%% select SVD filtering Noise
hhh=10;
load([IQfiles(hhh).folder filesep IQfiles(hhh).name],'IQ');
bulles = SVDfilter(IQ,ULM.SVD_cutoff);
bulles = filter(but_b,but_a,bulles,[],3);
bulles(~isfinite(bulles))=0;

figure(1)
dB = 20*log10(abs(bulles)); dB = dB-max(dB(:));
imagesc(lx,lz,dB(:,:,20),[-30 0]),colormap gray
colorbar,axis image,clear temp

BullesSNR = abs(bulles(:,:,10));
LocalMax = imregionalmax(BullesSNR);ValMax = sort(BullesSNR(LocalMax),'descend');
SNRmean = 20*log10(mean(ValMax(1:10))/mean(BullesSNR,'all'));clear BullesSNR LocalMax ValMax

%% Load and localize data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--- ULM PROCESSING --- \n\n')
clear Track_tot Track_tot_interp ProcessingTime bulles IQ dB
t1=tic;
parfor hhh = 1:min(Nbuffers,999) % can be run with for loop
    fprintf('Processing bloc %d/%d\n',hhh,Nbuffers);
    tmp = load([IQfiles(hhh).folder filesep IQfiles(hhh).name],'IQ');
    IQ_filt = SVDfilter(tmp.IQ,ULM.SVD_cutoff);tmp = [];
    IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0;

    % Data will be written in a '.mat' file to avoid RAM overdose
    [~,~] = PALA_multiULM(IQ_filt,listAlgo,ULM,PData,'savingfilename',[trackspath 'Tracks' num2str(hhh,'%.3d') '.mat']);
%         [Track_raw,Track_interp,ProTime] = PALA_multiULM(IQ_filt,listAlgo,ULM,PData);
%         save([trackspath 'Tracks' num2str(hhh,'%.3d')],'Track_raw','Track_interp','ProTime','ULM','UF','PData','-v6')
end
t2=toc(t1);
fprintf('ULM done in %d hours %.1f minutes. (for all localization algorithms) \n', floor(t2/60/60), rem(t2/60,60));

%% Create MatOuts     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each algorithm, create the MatOut density with interpolated tracks for visual analysis, and with non interpolated tracks for aliasing index calculation.
fprintf('--- CREATING MATOUTS --- \n\n')
MatOutSat = [];NbrOfLoc = zeros(Nalgo,1);ProcessingTime = zeros(Nbuffers,Nalgo);
MatOut = cell(Nalgo,1);MatOut(:)={0};MatOutNoInterp = MatOut;MatOut_vel = MatOut;

hwait = waitbar(0,'Bluiding intensity renderings','Name','Building matouts');
for hhh=1:min(Nbuffers,999) % Generate MatOut density matrix
    load([trackspath 'Tracks' num2str(hhh,'%.3d')],'Track_raw','Track_interp','ProTime')
    waitbar(hhh/Nbuffers,hwait);
    aa = -PData(1).Origin([3 1])+[1 1]*1;  % get origin
    bb = [1./PData(1).PDelta([3 1])*ULM.res];  % fix the size of pixel
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
save([savingpath 'MatOut_multi'],'MatOut','MatOut_vel','MatOutSat','NbrOfLoc','ULM','listAlgo','Nalgo','PData','UF');
save([savingpath 'MatOut_multi_nointerp'],'MatOutNoInterp','MatOutSat','NbrOfLoc','ULM','listAlgo','PData','Nalgo','UF','ProcessingTime');
close(hwait);clear hwait

t_end = toc(t_start);
fprintf('PALA_VivoBrain.m performed in %d hours and %.1f minutes\n',floor(t_end/60/60),rem(t_end/60,60));
fprintf('Next, please run PALA_VivoBrain_fig.m\n')
run('PALA_VivoBrain_fig.m')
return
%% Display MatOut Intensity
figure(90),clf
a = tight_subplot(2,ceil(numel(listAlgo)/2));
for ialgo=1:numel(listAlgo)
    axes(a(ialgo))
    imagesc(MatOut{ialgo}.^(1/3))
    axis image, colormap hot,caxis([0 8])
    title(listAlgo{ialgo},'interpreter','none')
end
linkprop(a,{'CLim','XLim','YLim'});
% saveas(gcf,[savingpath '_multiMatOut_disp .png'])

%% Display MatOut Velocimetry
wv = 1540/UF.TxFreq*1e-3*1e-3; %[mm]
figure(91),clf
a = tight_subplot(2,ceil(numel(listAlgo)/2));
for ialgo=1:numel(listAlgo)
    axes(a(ialgo))
    imagesc(MatOut_vel{ialgo}.*wv) % in [mm/s]
    axis image, colormap jet
    title(listAlgo{ialgo},'interpreter','none')
end
colorbar;linkprop(a,{'CLim','XLim','YLim'});
% saveas(gcf,[savingpath '_multiMatVel_disp .png'])
