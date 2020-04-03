function MatTracking = ULM_Superloc2D(MatIn,ULM)
%% function MatTracking = ULM_Superloc2D(MatIn,ULM)
% This function gets the localization for each particle using interpolation and loops
% This function was created by Baptiste Heiles and Ivan Piñeiro on 07/04/17
%
% INOUTS : 
%       - MatIn is the sequence containing all the images
%       - ULM structure must contains all parameters for localization method
%           - NLocalMax : number of local maxima
%           - LocMethod : localisation mehtod {'WA','Interp','Radial','CurveFitting','NoLocalization'};
%           - InterpMethod : {'bicubic','lanczos3','spline'} (methods avaliable for imresize)
%           - numberOfParticles is an estimation of the number of particle per image
%           - fwhm is the fwhm of a bubble (usually 3 if pixel at \lambda, 5 if pixel at \lambda/2)
% OUPUTS : 
%       - MatTracking is the table that stores the paricles values and position in [pixel]
% MatInReduced is the input matrix within a zeros frame of FWHM/2
%
% Created by Heiles Baptiste, Chavignon Arthur
% Last modifications Chavignon, 18/11/19
%
% DATE 2020.03.30 - VERSION 1.0.0
% AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
% Code Available under Creative Commons Non-Commercial 4.0
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018 


%% Get input data
fwhmz = ULM.fwhm(2);
fwhmx = ULM.fwhm(1);
numberOfParticles = ULM.numberOfParticles; 

%% Iinit struct
if ~isfield(ULM,'LocMethod')
    ULM.LocMethod = 'Radial';
end

if ~isfield(ULM,'parameters')
    ULM.parameters = struct();
end

if strcmp(ULM.LocMethod,'Interp')
    if ~isfield(ULM.parameters,'InterpMethod')
        ULM.parameters.InterpMethod = 'spline';
    end
    if sum(strcmp(ULM.parameters.InterpMethod,{'bilinear','bicubic'}))
        warning('Faster but very ugly.')
    end
end

if ~isfield(ULM.parameters,'NLocalMax')
    if fwhmz==3
        ULM.parameters.NLocalMax = 2;
    else
        ULM.parameters.NLocalMax = 3;
    end
end

%% Fill local parameters
NLocalMax = ULM.parameters.NLocalMax;
LocMethod = ULM.LocMethod;

%% MatIn format
info = whos('MatIn');
typename = info.class;

%% Initializing variables
[height,width,numberOfFrames]=size(MatIn);
MatIn = abs(MatIn);

%% Creates smaller matrix to avoid boundaries
MatInReduced = zeros(height,width,numberOfFrames,typename);
MatInReduced(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1, 1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:) = MatIn(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1, 1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:);
[height,width,numberOfFrames] = size(MatInReduced);

%% Imregionalmax implemented in 2D
Mat2D = permute(MatInReduced, [1,3,2]); %so that all the frames are in columns
Mat2D = reshape(Mat2D,height*numberOfFrames,width);
mask2D = imregionalmax(Mat2D,8); clear Mat2D
mask = reshape(mask2D,height,numberOfFrames,width);clear mask2D
mask = permute(mask,[1,3,2]);
IntensityMatrix = MatInReduced.*mask; %Values of intensities at regional maxima

%% Only keep highest values of Intensity
[tempMatrix,~] = sort(reshape(IntensityMatrix,[],size(IntensityMatrix,3)),1,'descend');

%% Preparing the intensities and coordinates for further calculation of average, intensities etc...
IntensityFinal = IntensityMatrix - ones(size(IntensityMatrix)) .* reshape(tempMatrix(numberOfParticles+1,:),[1 1 numberOfFrames]);
clear tempMatrix
MaskFinal   = (mask.*IntensityFinal)>0; %Mask with all intensities of bubbles low resolved
MaskFinal(isnan(MaskFinal))=0;
MaskFinal   = (MaskFinal>0).*IntensityMatrix;
index_mask  = find(MaskFinal);
[index_mask_z,index_mask_x,index_numberOfFrames]=ind2sub([height, width, numberOfFrames], index_mask);
clear mask IntensityFinal MaskFinal IntensityMatrix
clear index_mask

%% Creating FWHM models
% Building a vector from -FWHM to FWHM, this vector will be used for the mask's shifting
vectfwhmz = -1*round(fwhmz/2):round(fwhmz/2); 
vectfwhmx = -1*round(fwhmx/2):round(fwhmx/2); 

%% Localization
averageXc = nan(1,size(index_mask_z,1),typename);
averageZc = nan(1,size(index_mask_z,1),typename);

for iscat=1:size(index_mask_z,1)
    IntensityRoi  = (abs(MatIn(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat))));
    
    switch LocMethod
        case 'Radial'
            [Zc,Xc,sigma] = LocRadialSym(IntensityRoi,fwhmz,fwhmx);
        case 'WA'
            [Zc,Xc,sigma] = LocWeightedAverage(IntensityRoi,vectfwhmz,vectfwhmx);
        case 'Interp'
            [Zc,Xc,sigma] = LocInterp(IntensityRoi,ULM.parameters.InterpMethod,vectfwhmz,vectfwhmx);
        case 'CurveFitting'
            [Zc,Xc,sigma] = curveFitting(IntensityRoi,vectfwhmz,vectfwhmx);
        case 'NoLocalization'
            [Zc,Xc,sigma] = NoLocalization(IntensityRoi);
        otherwise
            error('Wrong LocMethod selected')
    end
    
    averageZc(iscat) = Zc + index_mask_z(iscat);
    averageXc(iscat) = Xc + index_mask_x(iscat);
    
    % implement in case you need to switch between models if the bubble is too large. For interp schemes there is liimted use to do this
    if or(sigma<0,sigma>25)
%         averageZc(iscat)=nan;
%         averageXc(iscat)=nan;
%         continue
    end
    
    % safeguard if too many localmax in IntensityROI
    if nnz(imregionalmax(IntensityRoi))>NLocalMax 
        averageZc(iscat)=nan;
        averageXc(iscat)=nan;
        continue
    end
    
    % if outside the center pixel of IntensityRoi
    if or(abs(Zc)>fwhmz/2,abs(Xc)>fwhmx/2) 
        averageZc(iscat)=nan;
        averageXc(iscat)=nan;
        continue
    end
end
keepIndex = ~isnan(averageXc);

ind = sub2ind([height,width,numberOfFrames],index_mask_z(keepIndex),index_mask_x(keepIndex),index_numberOfFrames(keepIndex));
clear index_mask_z index_mask_x IntensityRoi

%% Creating the table which stores the high resolved bubbles coordinates and the density value
MatTracking = zeros(nnz(keepIndex),4,typename);

MatTracking(:,1) = MatInReduced(ind);
MatTracking(:,2) = averageZc(keepIndex);
MatTracking(:,3) = averageXc(keepIndex);
MatTracking(:,4) = index_numberOfFrames(keepIndex);
clear averageXc averageZc index_numberOfFrames MatInReduced

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additionnal Functions

function sigma = ComputeSigmaScat(Iin,Zc,Xc)
[Nx,Nz] = size(Iin);
Isub = Iin - mean(Iin(:));
[px,pz] = meshgrid(1:Nx,1:Nz);
zoffset = pz - Zc+(Nz)/2.0;%BH xoffset = px - xc;
xoffset = px - Xc+(Nx)/2.0;%BH yoffset = py - yc;
r2 = zoffset.*zoffset + xoffset.*xoffset;
sigma = sqrt(sum(sum(Isub.*r2))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width
end

function [Zc,Xc,sigma] = LocRadialSym(Iin,fwhm_z,fwhm_x)
%% function [Zc,Xc,sigma] = LocRadialSym(Iin,fwhm_z,fwhm_x)
[Xc,Zc,sigma] = localizeRadialSymmetry(Iin,fwhm_z,fwhm_x);

end

function [Zc,Xc,sigma] = NoLocalization(Iin)
%% function [Zc,Xc,sigma] = NoLocalization(Iin)
Xc = 0;
Zc = 0;

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
Zc = sum(sum(Iin.*vectfwhm_z',1),2)./sum(Iin(:));
Xc = sum(sum(Iin.*vectfwhm_x,1),2)./sum(Iin(:));

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)

Nz=size(Iin,1);Nx=size(Iin,2);
if strcmp(InterpMode,'spline')
    [X,Z] = meshgrid(1:Nx,1:Nz);
    [Xq,Zq] = meshgrid(linspace(1,Nx,Nx*10),linspace(1,Nz,Nz*10));%to avoid uneven shift
    In_interp = interp2(X,Z,Iin,Xq,Zq,InterpMode);
else
    In_interp = imresize(Iin,10,InterpMode);
end
[~,tt] = max(In_interp(:));
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




