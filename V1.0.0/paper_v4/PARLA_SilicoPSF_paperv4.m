%% ULM PARLA: Post Processing - errors and pairings algorithms for IN SILICO PSF
% Perfoms ULM on Verasonics simulation for PSF variations.
% IQ are loaded, then noised at different clutter level. For each algorithm, the center
% position of the scatterer is computed with different algorithms and compare each others.
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
PARLA_data_folder = 'F:\ArthurC\Data\Fight_Club\Data';
addpath(genpath(PARLA_addons_folder))
%% Select IQ and Media folder
workingdir = [PARLA_data_folder '\PARLA_data_InSilicoPSF'];
filename = '2D_mesh_0205_120154';
cd(workingdir)

myfilepath = [workingdir filesep filename];
mkdir([workingdir filesep 'img'])
resultFile = [workingdir filesep 'img' filesep filename];

listVar = {'P','PData','Trans','Media','UF','Resource','filetitle','ll','mesh_x','mesh_z'};
load([myfilepath '_sequence.mat'],'-mat',listVar{:})

suffixe = '1';
%% ULM parameters
NFrames = P.BlocSize*P.numBloc;
res = 10;
ULM = struct('numberOfParticles', 1,...  % Number of particle per frame. (30-100)
    'res',10,... % Resolution factor. Typically 10 for images at lambda/10.
    'fwhm',[3 3],... % Size of the mask for localization. (3x3 pour BF à lambda, 5x5 pour lambda/2). [fmwhz fmwhx fmwhz]
    'size',[PData.Size(1),PData.Size(2),NFrames],...
    'scale',[1 1 1],... %sacle [z x y]
    'numberOfFramesProcessed',NFrames);%number of processed frames
res = ULM.res;

lx = PData.Origin(1) + [0:PData.Size(2)-1]*PData.PDelta(1);
lz = PData.Origin(3) + [0:PData.Size(1)-1]*PData.PDelta(3);

ULM.algorithm = 'wa'; % will be relace later by listAlgo(ii)

listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};
listAlgoName2 = {'No Shift','Weighted\newline Average',' Interp\newline Cubic',' Interp\newline Lanczos',' Interp\newline Spline','Gaussian\newline Fitting','Radial'};

Npoints = numel(ll)^2; % total number of simulated IQ (ie. bubble's position)
Nalgo = numel(listAlgo);

%% Simulated Noise paramters
% created a clutter-like noise.
NoiseParam.Power        = -2;   % [dBW]
NoiseParam.Impedance    = .2;   % [ohms]
NoiseParam.SigmaGauss   = 1.5;  % gaussian filtering
NoiseParam.clutterdB    = -30;  % clutter level in dB (will be changed later)
NoiseParam.amplCullerdB = 10;   % dB amplitude of clutter

% The code processes the simulated data for various SNR listed above :
listdB = [-60 -40 -30 -25 -20 -15 -10];  %dB

%% Display noised IQ example
temp = load([myfilepath '_IQ001.mat'],'IQ');
IQ_speckle = AddNoiseInIQ(abs(temp.IQ),NoiseParam);

figure(1)
imagesc(abs(IQ_speckle(:,:,2)))

%% Explaining nosie function
if 1
    IQ = abs(temp.IQ);
    IQ_speckle = AddNoiseInIQ(abs(temp.IQ),NoiseParam);
    
    ii=200;
    n_out = 10;n_in = size(IQ,1);
    IQ_0 = zeros(n_in+2*n_out,n_in+2*n_out);
    IQ_0(n_out + (0:n_in-1),n_out+(0:n_in-1)) = abs(IQ(:,:,ii));
    
    IQ_1 = reshape(wgn(numel(IQ_0),1,NoiseParam.Power,NoiseParam.Impedance),size(IQ_0,1),size(IQ_0,2),[]);
    IQ_1 = IQ_1*max(abs(IQ_0(:)))*10^((NoiseParam.amplCullerdB+NoiseParam.clutterdB)/20);
    IQ_2  = IQ_1+ max(abs(IQ_0(:)))*10^(NoiseParam.clutterdB/20);
    IQ_3 = imgaussfilt(IQ_2,NoiseParam.SigmaGauss);
    IQ_4 = IQ_3+IQ_0;
    IQ_dB = 20*log10(abs(IQ_4));IQ_dB = IQ_dB -max(IQ_dB (:));
    figure(101);clf;
    subplot 241
    imagesc(abs(IQ_0)),axis image,axis off,title('raw');colorbar
    subplot 242
    imagesc(abs(IQ_1)),axis image,axis off,title('noise normalized');colorbar
    subplot 243
    imagesc(abs(IQ_2)),axis image,axis off,title('noise +clutter');colorbar
    subplot 244
    imagesc(abs(IQ_3)),axis image,axis off,title('noise gauss');colorbar
    subplot 245
    imagesc(abs(IQ_4)),axis image,axis off,title('Noised PSF');colorbar
    subplot 246
    imagesc(IQ_dB,[-60 0]),colorbar;axis image,axis off
    subplot 247
    hold off
    plot(IQ_0(:,4))
    hold on
    plot(IQ_1(:,4))
    plot(IQ_2(:,4))
    plot(IQ_3(:,4))
    plot(IQ_4(:,4),'k:','LineWidth',2)
    subplot 248
    plot(IQ_dB)
    ylim([-50 0])
    for ii=1:size(IQ_4,3),SNR_measure(ii) = 20*log10(max(reshape(IQ_4(:,:,ii),[],1)) / mean(reshape(IQ_4(8:end,:,ii),[],1)));end
    title(mean(SNR_measure))
end

%% Display IQnoise and scatter's position
dB = 20*log10(abs(IQ_speckle)); dB = dB-max(dB(:));

figure(1)
for ii=1:10:100%size(IQ_speckle,3)
    imagesc(lx,lz,dB(:,:,ii),[-30 0]),colormap gray;
    hold on
    colorbar,axis image
    plot(Media.ListPos(ii,1),Media.ListPos(ii,3),'r.')
    drawnow;pause(0.01)
end

%% Load and localize data

for idB = 1:numel(listdB)
    %% Simulation and loczalization for each clutterLevel
    NoiseParam.clutterdB = listdB(idB)
    MatLocFull = []; %[pos algo frame rep]
    load([myfilepath '_IQ001.mat'],'IQ');
    
    % Number of repetition of noise occurence (for each position, Nt simulation of noise are performs to randomize effects.
    Nt = 2;
    for irp = 1:Nt
        %% Perform localization using a dedicated function FC_Superloc_multi2_mesh
        fprintf(['\n ' num2str(irp) ': '])
        IQ_speckle = AddNoiseInIQ(abs(IQ),NoiseParam);
        MatTracking = FC_Superloc_multi2_mesh(IQ_speckle,listAlgo,ULM,PData);
        MatLocFull = cat(4,MatLocFull,MatTracking); % store data Dimension [position / list algo / frames / repetition]
    end
    % Store and save raw localizations
    save([myfilepath '_LocalMesh'  num2str(abs(NoiseParam.clutterdB)) 'dB'],'MatLocFull','Media','Nt','ULM','P','listAlgo','NoiseParam','ll','listAlgoName2');

    %% Reshape and calulate errors
    MatLocFull_w =    MatLocFull; % create a copy    
    MatPosSim =     permute(Media.ListPos(:,[3 1]),[2 3 1]); % theorical position's of scatterers
    MatLocFull_w =    MatLocFull_w(:,:,1:size(MatPosSim,3),:); % remove out of list positions
    ErrorFull =     MatLocFull_w-MatPosSim; % compute localization errors
    
    % Data are reshape to get a user-friendly organization.
    % for each of the Npoints [21x21 grid] position, we have a 3 random shift sub position (ie [ix/21 + rand/21, iz/21 + rand/21).
    % The calcultation is computed fo each Algo (actually 7), and for Nt occasion of noise.
    % The error is sorted for anaylsing result in the 21x21 grid -> [pos algo GridPos repetitions]
    
    ErrorFull =     reshape(ErrorFull,2,Nalgo,Npoints,[],Nt); % Dim [pos, algo, grid point, random shift, noise repetition]
    ErrorFull =     permute(ErrorFull,[1 2 3 5 4]); % Dim [pos, algo, grid point, noise repetition, random shift]
    ErrorFull =     reshape(ErrorFull,2,Nalgo,Npoints,[]); % Dim [pos, algo, grid point, total number of repetition]
    
    %% Remove static error
    % It remains a small localization error which is specific for each algorihm. It can be
    % callibrated by measuring the error with a very low clutter level. Here we estimated
    % this erros a -60dB. The value is store for each algorihm at -60dB, and values are
    % corrected later.
    ListMean(:,:,:) = mean(ErrorFull(1:2,:,:,:),4);
    ListMean(3,:,:) = mean(sqrt(sum(ErrorFull(1:2,:,:,:).^2,1)),4);
    
    if abs(NoiseParam.clutterdB)==60
        StaticError = mean(ListMean,3); % measure the specific localization for each algo.
    else
        % load the specific localization at -60dB
        load([myfilepath '_LocalMesh'  num2str(abs(-60)) 'dB'],'StaticError');
    end
    % Correct localization errors with this specific static error.
    ListMean = ListMean-StaticError;    
    ListMean(3,:,:) = mean(sqrt(sum((ErrorFull(1:2,:,:,:)-StaticError(1:2,:,1,1)).^2,1)),4);
    
    ListVar = var(ErrorFull(1:2,:,:,:),0,4);
    ListVar(3,:,:) = var(sqrt(sum(ErrorFull(1:2,:,:,:).^2,1)),0,4);
    
    save([myfilepath '_LocalMesh'  num2str(abs(NoiseParam.clutterdB)) 'dB'],'-append','StaticError','ListMean');
    
    %% Display targets and localizations
    figure(13),clf;tightax=tight_subplot(2,4,0.03,0.03,0.01);
    for ialgo=1:Nalgo
        axes(tightax(ialgo)),hold on
        plot(Media.ListPos(:,1),Media.ListPos(:,3),'k.'),axis image;xlabel('\lambda');ylabel('\lambda')
        plot(squeeze(MatLocFull_w(2,ialgo,:,1))-StaticError(2,ialgo),squeeze(MatLocFull_w(1,ialgo,:,1))-StaticError(1,ialgo),'b+')
        title(listAlgo(ialgo),'Interpreter','none')
    end
%     saveas(gcf,[resultFile '_TargetNShot' num2str(NoiseParam.clutterdB) 'dB.png'])
    
    %% Display Errors maps
    Caxis.x = [0 1]*.5;
    Caxis.y = [0 1]*.5;
    Caxis.RMSE = [0 1]*.5;
    
    f101 = figure(101);f101.Position = [510 521 1170 457];clf
    colormap jet
    
    aa = tight_subplot(3,7,[.05 .01],[.05 .09],[.04 .05]);
    
    ErrorMaps = abs(ListMean);
    ErrorMaps = reshape(ErrorMaps,[3 Nalgo numel(ll) numel(ll)]);
    ErrorMaps = permute(ErrorMaps,[3 4 1 2]);
    
    for ialgo = 1:Nalgo
        axes(aa(ialgo))
        imagesc(llx,llz,ErrorMaps(:,:,2,ialgo)),axis image
        caxis(Caxis.x)
        title(listAlgoName2(ialgo))
        axis off
        if ialgo==1
            ylabel('Lateral Error');
            axis on;
        end
    end
    clb = colorbar;clb.Position(1)=aa(7).Position(1)+aa(7).Position(3);clb.Position([2 4]) = aa(7).Position([2 4]);
    clb.Label.String='Absolute error';
    
    for ialgo = 1:Nalgo
        axes(aa(ialgo+Nalgo))
        imagesc(llx,llz,ErrorMaps(:,:,1,ialgo)),axis image
        caxis(Caxis.y)
        axis off
        if ialgo==1,ylabel('Axial Error');axis on;end
    end
    ylabel('Z')
    cbr_z = colorbar;cbr_z.Position(1) = aa(7+Nalgo).Position(1) + aa(7+Nalgo).Position(3);cbr_z.Position([2 4]) = aa(7+Nalgo).Position([2 4]);
    cbr_z.Label.String='Absolute error';
    
    for ialgo = 1:Nalgo
        axes(aa(ialgo+Nalgo*2))
        imagesc(llx,llz,ErrorMaps(:,:,3,ialgo)),axis image
        caxis(Caxis.RMSE)
        axis off
        if ialgo==1,ylabel('RMSE');axis on;end
    end
    cbr_r = colorbar;cbr_r.Position(1) = aa(7+Nalgo*2).Position(1) + aa(7+Nalgo*2).Position(3);cbr_r.Position([2 4]) = aa(7+Nalgo*2).Position([2 4]);
    cbr_r.Label.String='Absolute error';
%     saveas(gcf,[resultFile '_MeanErrorMap' num2str(NoiseParam.clutterdB) 'dB.png'])

end

