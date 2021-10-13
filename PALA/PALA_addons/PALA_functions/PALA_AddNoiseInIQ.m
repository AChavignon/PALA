function IQ_speckle = PALA_AddNoiseInIQ(IQ,NoiseParam)
%% function IQ_speckle = PALA_AddNoiseInIQ(IQ,NoiseParam)
% Takes raw IQ in input and returns a noised IQ simulating a clutter noise.
% requires Communications Toolbox.
%
% Created by Arthur Chavignon 25/02/2020
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018
%
% Noise Parameters example
%       - NoiseParam.Power        = -2; %[dBW]
%       - NoiseParam.Impedance    = .2; %[ohms]
%       - NoiseParam.SigmaGauss   = 1.5; % Gaussian filtering
%       - NoiseParam.clutterdB    = -20; % clutter level in dB
%       - NoiseParam.amplCullerdB = 10; % dB amplitude of clutter

IQ = abs(IQ);
IQ_speckle = IQ+imgaussfilt(max(IQ(:))*10^(NoiseParam.clutterdB/20)+reshape(wgn(numel(IQ),1,NoiseParam.Power,NoiseParam.Impedance),size(IQ,1),size(IQ,2),[])*max(IQ(:))*10^((NoiseParam.amplCullerdB+NoiseParam.clutterdB)/20),NoiseParam.SigmaGauss);

%% Details:
if 0
    ii=max(200,size(IQ,3));
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
    subplot 241,imagesc(abs(IQ_0)),axis image,axis off,title('raw');colorbar
    subplot 242,imagesc(abs(IQ_1)),axis image,axis off,title('noise normalized');colorbar
    subplot 243,imagesc(abs(IQ_2)),axis image,axis off,title('noise +clutter');colorbar
    subplot 244,imagesc(abs(IQ_3)),axis image,axis off,title('noise gauss');colorbar
    subplot 245,imagesc(abs(IQ_4)),axis image,axis off,title('Noised PSF');colorbar
    subplot 246,imagesc(IQ_dB,[-60 0]),colorbar;axis image,axis off
    subplot 247,hold off
    plot(IQ_0(:,4))
    hold on
    plot(IQ_1(:,4));plot(IQ_2(:,4));plot(IQ_3(:,4));plot(IQ_4(:,4),'k:','LineWidth',2)
    subplot 248,plot(IQ_dB)
    ylim([-50 0])
    for ii=1:size(IQ_4,3),SNR_measure(ii) = 20*log10(max(reshape(IQ_4(:,:,ii),[],1)) / mean(reshape(IQ_4(8:end,:,ii),[],1)));end
    title(mean(SNR_measure))
end

end
