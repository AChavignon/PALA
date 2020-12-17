%% PALA_SetUpPaths.m : define path for data and addons
% Before running scripts, please fill data and addons paths.
%
% Created by Arthur Chavignon 17/12/2020
%
% DATE 2020.12.17 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

clear all;close('all')

% DEFINE THE ADDONS DIRECTORY ON YOUR COMPUTER

PALA_addons_folder = 'F:\ArthurC\CHAPRO_local\_SIMULATIONS\FightClub_PostPro\PALA_addons'; % location of the addons folder

% DEFINE THE DATA DIRECTORY ON YOUR COMPUTER

PALA_data_folder = 'F:\ArthurC\Data\Fight_Club\PALA_Data';

addpath(genpath(PALA_addons_folder))

%% Checking licenses and features
ToolBoxRequires = {'Communications','Bioinformatics','Image Processing','Curve Fitting','Signal Processing','Statistics and Machine Learning','Parallel Computing','Computer Vision Toolbox'};
err = 0;
for featureName=ToolBoxRequires
   IsInstalledToolbox = contains(struct2array(ver), featureName{1});
   if ~IsInstalledToolbox, warning([featureName{1} ' is missing']),err=1;end
end
if err,error('Toolbox are missing.');end;clear ToolBoxRequires featureName IsInstalledToolbox err