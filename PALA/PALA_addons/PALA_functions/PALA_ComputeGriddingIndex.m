function [Gridding_index,varargout] = PALA_ComputeGriddingIndex(MatOutNoInterp,i_harmonic)
%% function Gridding_index = PALA_ComputeGriddingIndex(MatOutNoInterp,i_harmonic)
% This function computes the Gridding index of non interpolated rendering of ULM.
% In the Fourier Space, it compares the intensity of peaks to the baseline. Peaks
% appears when the initial beamforming grid appears in the MatOut (ie. when localization
% are kept in the center of the beamforming grid).
%
% INPUTS:
%       - MatOutNoInterp : cell of MatOut with non-interpolated tracks. (ie. only localization)
%       - i_harmonic : list of harmonics used for the calculation (default [1 2])
% OUTPUTS :
%       - Gridding_index : gridding index
%       - varargout{1}=F_0: fondamentala frequency of the gridding effect
%
% Created by Arthur Chavignon 9/12/2019
%
% DATE 2020.08.10 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

if ~exist('i_harmonic')
    i_harmonic = [1 2]';
end

i_harmonic = reshape(i_harmonic,[],1);
Nharm = numel(i_harmonic);
%% Find first harmonic in x and z directions
res = 10; % size of SRpixel compare to original one
F_0 = round(size(MatOutNoInterp{1})/res); % frequency of aliasing (related to the size of initial pixel of the image)

disp(['Fundamental frequencies: ' num2str(F_0)])
if nargout>1
    varargout{1}=F_0;
end

%% Frequencies to sum
Freq_i_z = F_0(1)*i_harmonic + reshape([-1:1],1,[]);
Freq_i_x = F_0(2)*i_harmonic + reshape([-1:1],1,[]);

%% Frequencies for the baseline
Freq_base_z = F_0(1)*i_harmonic + [-5:-3,3:5];
Freq_base_x = F_0(2)*i_harmonic + [-5:-3,3:5];

%% Compute of each matout
TF2_z_list=[]; TF2_x_list=[];
Gridding_index = zeros(1,numel(MatOutNoInterp));
Gridding_x = Gridding_index;Gridding_z = Gridding_index;

for ii=1:numel(MatOutNoInterp)
    %% 2D Fourier Transform of the MatOut non interpolated
    TF2 = (fft2(MatOutNoInterp{ii}));
    
    % projection along x and z direction, and normalization
    TF2_z = smooth(sum(abs(TF2).^2,2)/sum(sum(abs(TF2).^2,1)),1);
    TF2_x = smooth(sum(abs(TF2).^2,1)/sum(sum(abs(TF2).^2,1)),1);
    
    %     TF2_z = TF2_z(1:(size(TF2_z,1)-1)/2) + flip( TF2_z((size(TF2_z,1)-1)/2+2:end) );
    %     TF2_x = TF2_x(1:(size(TF2_x,1)-1)/2) + flip( TF2_x((size(TF2_x,1)-1)/2+2:end) );
    
    peak_z = max(reshape(TF2_z(Freq_i_z),Nharm,[]) ,[],2);
    baseline_z = mean(reshape(TF2_z(Freq_base_z),Nharm,[]),2);
    Gridding_z(ii) = 20*log10(sum(peak_z./baseline_z));
    
    peak_x = max(reshape(TF2_x(Freq_i_x),Nharm,[]),[],2);
    baseline_x = mean(reshape(TF2_x(Freq_base_x),Nharm,[]),2);
    Gridding_x(ii) = 20*log10(sum(peak_x./peak_x));
    
    Gridding_index(ii) = 20*log10(1/2*sum(peak_x./baseline_x  + peak_z./baseline_z));
    
    TF2_z_list(ii,:) = TF2_z;
    TF2_x_list(ii,:) = TF2_x;
    
end

end
