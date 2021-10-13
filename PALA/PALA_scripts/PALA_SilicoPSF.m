%% PALA_SilicoPSF.m : Post Processing - errors and pairings algorithms for IN SILICO PSF
% Performs ULM on Verasonics simulation for PSF variations.
% IQ are loaded, then noised at different clutter level. For each algorithm, the center
% position of the scatterer is computed with different algorithms and compare each others.
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
fprintf('Running PALA_SilicoPSF.m\n')
workingdir = [PALA_data_folder '\PALA_data_InSilicoPSF'];cd(workingdir)
filename = 'PALA_InSilicoPSF';
myfilepath = [workingdir filesep filename];

listVar = {'P','PData','Trans','Media','UF','Resource','filetitle','ll','mesh_x','mesh_z'};
load([myfilepath '_sequence.mat'],'-mat',listVar{:})
%% ULM parameters
NFrames = P.BlocSize*P.numBloc;
ULM = struct('numberOfParticles', 1,...  % Number of particles per frame. (30-100)
    'res',10,... % Resolution factor. Typically, 10 for images at lambda/10.
    'fwhm',[3 3],... % Size of the mask for localization. (3x3 for pixel at lambda, 5x5 at lambda/2). [fmwhz fmwhx]
    'size',[PData.Size(1),PData.Size(2),NFrames],...
    'scale',[1 1 1],... % Scale [z x y]
    'algorithm','wa',... % Will be replaced later by listAlgo(ii)
    'numberOfFramesProcessed',NFrames); % Number of processed frames
res = ULM.res;

lx = PData.Origin(1) + [0:PData.Size(2)-1]*PData.PDelta(1);
lz = PData.Origin(3) + [0:PData.Size(1)-1]*PData.PDelta(3);

listAlgo = {'no_shift','wa','interp_cubic','interp_lanczos','interp_spline','gaussian_fit','radial'};

Npoints = numel(ll)^2; % Total number of simulated IQ (ie. bubble's position)
Nalgo = numel(listAlgo);

%% Simulated Noise parameters'
% created a clutter-like noise.
NoiseParam.Power        = -2;   % [dBW]
NoiseParam.Impedance    = .2;   % [ohms]
NoiseParam.SigmaGauss   = 1.5;  % Gaussian filtering
NoiseParam.clutterdB    = -30;  % Clutter level in dB (will be changed later)
NoiseParam.amplCullerdB = 10;   % dB amplitude of clutter

% The code processes the simulated data for various SNR listed above:
listdB = [-60 -40 -30 -25 -20 -15 -10];  %dB
%% Display noised IQ example
temp = load([myfilepath '_IQ001.mat'],'IQ');
IQ_speckle = PALA_AddNoiseInIQ(abs(temp.IQ),NoiseParam);
figure(1);imagesc(abs(IQ_speckle(:,:,2)))

%% Display IQnoise and scatter's position
dB = 20*log10(abs(IQ_speckle)); dB = dB-max(dB(:));
figure(1)
for ii=1:10:100 %size(IQ_speckle,3)
    imagesc(lx,lz,dB(:,:,ii),[-30 0]),colormap gray;hold on;colorbar,axis image
    plot(Media.ListPos(ii,1),Media.ListPos(ii,3),'r.')
    drawnow;pause(0.01)
end

%% Load and localize data
t_start=tic;
for idB = 1:numel(listdB)
    % Simulation and loczalization for each clutterLevel
    disp(['Simulation ' num2str(idB) '/' num2str(numel(listdB)) ': SNR = ' num2str(abs(listdB(idB))) 'dB'])
    clear ListMean MatLocFull_w MatPosSim ErrorFull ListVar
    NoiseParam.clutterdB = listdB(idB);
    MatLocFull = []; %[pos algo frame rep]
    load([myfilepath '_IQ001.mat'],'IQ');

    % Number of repetition of noise occurrence (for each position, Nt simulation of noise are performs to randomize effects.
    Nt = 2;
    for irp = 1:Nt
        % Perform localization using a dedicated function PALA_multiULM_mesh
        fprintf(['\n Iteration ' num2str(irp) ': '])
        IQ_speckle = PALA_AddNoiseInIQ(abs(IQ),NoiseParam);
        MatTracking = PALA_multiULM_mesh(IQ_speckle,listAlgo,ULM,PData);
        MatLocFull = cat(4,MatLocFull,MatTracking); % store data Dimension [position / list algo / frames / repetition]
    end
    % Store and save raw localizations
    save([myfilepath '_LocalMesh'  num2str(abs(NoiseParam.clutterdB)) 'dB'],'MatLocFull','Media','Nt','ULM','P','listAlgo','NoiseParam','Nalgo','ll');

    %% Reshape and compute errors
    MatLocFull_w =  MatLocFull; % create a copy
    MatPosSim =     permute(Media.ListPos(:,[3 1]),[2 3 1]); % true position's of scatterers
    MatLocFull_w =  MatLocFull_w(:,:,1:size(MatPosSim,3),:); % remove out of list positions
    ErrorFull =     MatLocFull_w-MatPosSim; % compute localization errors

    % Data are reshape to get a user-friendly organization.
    % for each of the Npoints [21x21 grid] position, we have a 3 random shift sub position (ie [ix/21 + rand/21, iz/21 + rand/21).
    % The calculation is computed for each algorithms (actually 7), and for Nt occasion of noise.
    % The error is sorted for analyzing result in the 21x21 grid -> [pos algo GridPos repetitions]
    ErrorFull =     reshape(ErrorFull,2,Nalgo,Npoints,[],Nt); % Dim [pos, algo, grid point, random shift, noise repetition]
    ErrorFull =     permute(ErrorFull,[1 2 3 5 4]); % Dim [pos, algo, grid point, noise repetition, random shift]
    ErrorFull =     reshape(ErrorFull,2,Nalgo,Npoints,[]); % Dim [pos, algo, grid point, total number of repetition]

    %% Remove static error
    % It remains a small localization error which is specific for each algorithm. It can be
    % calibrated by measuring the error with a very low clutter level. Here we estimated
    % this error at -60dB. The value is store for each algorithm at -60dB, and values are
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

    %% Display Errors' maps
    Caxis.x = [0 1]*.5;Caxis.y = [0 1]*.5;Caxis.RMSE = [0 1]*.5;
    f101 = figure(101);f101.Position = [510 521 1170 457];clf;colormap jet
    aa = tight_subplot(3,7,[.05 .01],[.05 .09],[.04 .05]);
    ErrorMaps = abs(ListMean);ErrorMaps = reshape(ErrorMaps,[3 Nalgo numel(ll) numel(ll)]);
    ErrorMaps = permute(ErrorMaps,[3 4 1 2]);

    for ialgo = 1:Nalgo
        axes(aa(ialgo));imagesc(ll,ll,ErrorMaps(:,:,2,ialgo)),axis image,caxis(Caxis.x),axis off
        title(listAlgo(ialgo),'interpreter','none')
        if ialgo==1,ylabel('Lateral Error');axis on;end
    end
    clb = colorbar;clb.Position(1)=aa(7).Position(1)+aa(7).Position(3);clb.Position([2 4]) = aa(7).Position([2 4]);
    clb.Label.String='Absolute error';

    for ialgo = 1:Nalgo
        axes(aa(ialgo+Nalgo));imagesc(ll,ll,ErrorMaps(:,:,1,ialgo)),axis image,caxis(Caxis.y);axis off
        if ialgo==1,ylabel('Axial Error');axis on;end
    end
    ylabel('Z')
    cbr_z = colorbar;cbr_z.Position(1) = aa(7+Nalgo).Position(1) + aa(7+Nalgo).Position(3);cbr_z.Position([2 4]) = aa(7+Nalgo).Position([2 4]);
    cbr_z.Label.String='Absolute error';

    for ialgo = 1:Nalgo
        axes(aa(ialgo+Nalgo*2));imagesc(ll,ll,ErrorMaps(:,:,3,ialgo)),axis image,caxis(Caxis.RMSE);axis off
        if ialgo==1,ylabel('RMSE');axis on;end
    end
    cbr_r = colorbar;cbr_r.Position(1) = aa(7+Nalgo*2).Position(1) + aa(7+Nalgo*2).Position(3);cbr_r.Position([2 4]) = aa(7+Nalgo*2).Position([2 4]);
    cbr_r.Label.String='Absolute error';drawnow
%     saveas(gcf,[myfilepath '_MeanErrorMap' num2str(NoiseParam.clutterdB) 'dB.png'])

end
t_end = toc(t_start);
fprintf('PALA_SilicoPSF.m performed in %d hours and %.1f minutes\n',floor(t_end/60/60),rem(t_end/60,60));
fprintf('Next, please run PALA_SilicoPSF_fig.m\n')
run('PALA_SilicoPSF_fig.m')

