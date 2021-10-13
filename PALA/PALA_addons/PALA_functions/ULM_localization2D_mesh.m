function MatTracking = ULM_localization2D_mesh(MatIn,ULM)
%% function MatTracking = ULM_localization2D_mesh(MatIn,ULM)
% This function performs the detection, selection and sub-pixel localization of bubbles on a list of input images (MatIn).
%
% MatIn is the sequence containing all the images
% fwhm is the fwhm of a bubble (usually 2 for 9MHz)
% numberOfParticles is an estimation of the number of particle per image
% MatTracking is the table that stores the paricles values and position
% MatOut is the Super-resolved image mask is the imregionalmax result. It is a logical matrix
% MatInReduced is the input matrix within a zeros frame of FWHM/2
%
% INPUTS:
%       - MatIn is the sequence containing all the images
%       - ULM structure must contains all parameters for localization method
%           - NLocalMax : number of local maxima
%           - LocMethod : localisation mehtod {'WA','Interp','Radial','CurveFitting','NoLocalization'};
%           - InterpMethod : {'bicubic','lanczos3','spline'} (methods avaliable for imresize)
%           - fwhm is the fwhm of a bubble (usually 3 if pixel at \lambda, 5 if pixel at \lambda/2)
% OUPUT:
%       - MatTracking is the table that stores the paricles values and position in [pixel]
%
% This function was created by Baptiste Heiles 07/04/17, last modifications Arthur Chavignon, 02/10/19
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

%% Iinit struct
if ~isfield(ULM,'LocMethod'),ULM.LocMethod = 'Radial';end
if ~isfield(ULM,'parameters'),ULM.parameters = struct();end

if strcmp(ULM.LocMethod,'Interp')
    if ~isfield(ULM.parameters,'InterpMethod')
        ULM.parameters.InterpMethod = 'spline';
    end
    if sum(strcmp(ULM.parameters.InterpMethod,{'bilinear','bicubic'}))
        warning('Faster but very ugly.')
    end
end

%% Fill parameters
LocMethod = ULM.LocMethod;
fwhmz = ULM.fwhm(2); fwhmx = ULM.fwhm(1);

%% Initializing variables
[height,width,numberOfFrames]=size(MatIn);
MatIn = abs(MatIn);

MatInReduced = zeros(height,width,numberOfFrames);
redx = 3:width-2;
redz = 3:height-2;
MatInReduced(redx, redz,:) = MatIn(redx,redz,:);

MatIn2D = reshape(MatInReduced,height*width,numberOfFrames);
[~,index_mask] =max(MatIn2D,[],1);
[index_mask_z,index_mask_x]=ind2sub([height, width], index_mask);

%% Creating FWHM models
% Building a vector from -FWHM to FWHM, this vector will be used for the mask's shifting
vectfwhmz = -1*round(fwhmz/2):round(fwhmz/2);
vectfwhmx = -1*round(fwhmx/2):round(fwhmx/2);

%% Localization
averageXc = nan(1,size(index_mask_z,1));
averageZc = nan(1,size(index_mask_z,1));

for iscat=1:size(MatIn,3)
    IntensityRoi  = MatIn(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,iscat);
    
    switch LocMethod
        case 'Radial'
            [Zc,Xc,~] = LocRadial(IntensityRoi,fwhmz,fwhmx);
        case 'WA'
            [Zc,Xc,~] = LocWeightedAverage(IntensityRoi,vectfwhmz,vectfwhmx);
        case 'Interp'
            [Zc,Xc,~] = LocInterp(IntensityRoi,ULM.parameters.InterpMethod,vectfwhmz,vectfwhmx);
        case 'CurveFitting'
            [Zc,Xc,~] = curveFitting(IntensityRoi,vectfwhmz,vectfwhmx);
        case 'NoLocalization'
            [Zc,Xc,~] = NoLocalization(IntensityRoi);
        otherwise
            error('Wrong LocMethod selected')
    end
    
    averageZc(iscat) = Zc + index_mask_z(iscat);
    averageXc(iscat) = Xc + index_mask_x(iscat);
end
clear index_mask_z index_mask_x

%% Creating the table which stores the high resolved bubbles coordinates and the density value
MatTracking = cat(2,averageZc(:),averageXc(:));
clear averageXc averageZc index_numberOfFrames MatInReduced

end

%% Additional Functions

function sigma = ComputeSigmaScat(Iin,Zc,Xc)
[Nx,Nz] = size(Iin);
Isub = Iin - mean(Iin(:));
[px,pz] = meshgrid(1:Nx,1:Nz);
zoffset = pz - Zc+(Nz)/2.0;%BH xoffset = px - xc;
xoffset = px - Xc+(Nx)/2.0;%BH yoffset = py - yc;
r2 = zoffset.*zoffset + xoffset.*xoffset;
sigma = sqrt(sum(sum(Isub.*r2))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width
end

function [Zc,Xc,sigma] = LocRadial(Iin,fwhm_z,fwhm_x)
%% function [Zc,Xc,sigma] = LocRadial(Iin,fwhm_z,fwhm_x)
[Zc,Xc] = localizeRadialSymmetry(Iin,fwhm_z,fwhm_x);
sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = NoLocalization(Iin)
%% function [Zc,Xc,sigma] = NoLocalization(Iin)
Xc = 0;
Zc = 0;

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
%%function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
Zc = sum(sum(Iin.*vectfwhm_z',1),2)./sum(Iin(:));
Xc = sum(sum(Iin.*vectfwhm_x,1),2)./sum(Iin(:));

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)
%%function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)

Nz=size(Iin,1);Nx=size(Iin,2);
if strcmp(InterpMode,'spline')
    [X,Z] = meshgrid(1:Nx,1:Nz);
    [Xq,Zq] = meshgrid(linspace(1,Nx,Nx*10),linspace(1,Nz,Nz*10));%to avoid uneven shift
    In_interp = interp2(X,Z,Iin,Xq,Zq,InterpMode);
else
    %     if strcmp(InterpMode,'lanczos3')
    %         In_interp = imresize(Iin,10,{@mylanczos2,7});
    %     else
    In_interp = imresize(Iin,10,InterpMode);
    %     end
end
[~,tt] = max(In_interp(:));
% [tt]=find(abs(In_interp)==max(abs(In_interp(:))));
% if size(tt,1)>1
%     tt=tt(1,1);%arbitrary solve this
% end
[iz,ix,~]=ind2sub(size(In_interp),tt);

Zc = vectfwhm_z(1)-1 + 1 + iz./10 -.5 +.05;
Xc = vectfwhm_x(1)-1 + 1 + ix./10 -.5 +.05;

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
[meshX,meshZ] = meshgrid(vectfwhm_x,vectfwhm_z);
meshIn = cat(3,meshX,meshZ);

sigGauss_z = vectfwhm_z(end)*0+1;
sigGauss_x = vectfwhm_x(end)*0+1;

myGaussFunc = @(x_pos,mesh_pos)( exp(-(mesh_pos(:,:,1)-x_pos(1)).^2./(2*sigGauss_z^2) - (mesh_pos(:,:,2)-x_pos(2)).^2./(2*sigGauss_x^2)));
OPTIONS = optimoptions('lsqcurvefit','StepTolerance',.01,'MaxIterations',5,'Display','off');

% Gaussian Fitting
x_out = lsqcurvefit(myGaussFunc,[0 0],meshIn,double(Iin./max(Iin(:))),[],[],OPTIONS);

Zc = x_out(2);
Xc = x_out(1);

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function f = mylanczos2(x)
f = (sin(pi*x) .* sin(pi*x/2) + eps) ./ ((pi^2 * x.^2 / 2) + eps);
f = f .* (abs(x) < 2);
end
