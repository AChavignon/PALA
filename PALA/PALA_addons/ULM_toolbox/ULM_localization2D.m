function MatTracking = ULM_localization2D(MatIn,ULM)
%% function MatTracking = ULM_localization2D(MatIn,ULM)
% This function performs the detection, selection and sub-pixel localization of bubbles on
% a list of input images (MatIn).
%
% - The detection step is performed with the imregionalmax function. It returns the list of local maxima.
% - The selection step consists of sorting intensities, in each frames, and keeping the highest maxima.
% - The localization steps consists of applying a sub-wavelength localization kernel
% (weighted average, interpolation, radial symmetry...) to a cropped image centered on a
% local maxima (a 5x5 or 3x3 large image). Localization kernel are discussed in the cited article.
%
% This function can be easily adapt by anyone who wants to try a new localization kernel,
% keeping all the framework unchanged.
%
% INPUTS:
%       - MatIn is the sequence containing all the images
%       - ULM structure must contains all parameters for localization method
%           - NLocalMax : number of local maxima
%           - LocMethod : localisation mehtod {'wa','interp','radial','curvefitting','nolocalization'};
%           - InterpMethod : {'bicubic','lanczos3','spline'} (methods avaliable for imresize)
%           - numberOfParticles is an estimation of the number of particles per image
%           - fwhm is the fwhm of a bubble (usually 3 if pixel at \lambda, 5 if pixel at \lambda/2)
% OUTPUT:
%       - MatTracking is the table that stores the particles values and position in [pixel]
% MatInReduced is the input matrix within a zeros frame of FWHM/2
%
% This function was created by Baptiste Heiles 07/04/17, last modifications Arthur Chavignon, 18/11/19
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Baptiste Heiles, Arthur Chavignon, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

%% Get input data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwhmz = ULM.fwhm(2);
fwhmx = ULM.fwhm(1);

% Building a vector from -FWHM to FWHM, this vector will be used for the mask's shifting
vectfwhmz = -1*round(fwhmz/2):round(fwhmz/2);
vectfwhmx = -1*round(fwhmx/2):round(fwhmx/2);

[height,width,numberOfFrames]=size(MatIn);% Get sizes of the Matrix, height denotes number of rows/depth of imaging, width denotes number of lines/width of imaging,
% numberOfFrames denotes the number of elements in the third dimension/number of Frames
MatIn = abs(MatIn);% Make sure you work with the intensity matrix
info = whos('MatIn');typename = info.class;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize structures
% if fields are missing, they will be set with default values.
if ~isfield(ULM,'LocMethod'),ULM.LocMethod = 'radial';end

if ~isfield(ULM,'parameters')% Create an empty structure for parameters hosting
    ULM.parameters = struct();
end

if strcmp(ULM.LocMethod,'interp')
    if ~isfield(ULM.parameters,'InterpMethod')
        ULM.parameters.InterpMethod = 'spline';
    end
    if sum(strcmp(ULM.parameters.InterpMethod,{'bilinear','bicubic'}))
        warning('Faster but pixelated, Weighted Average will be faster and smoother.')
    end
end

if ~isfield(ULM.parameters,'NLocalMax')
    if fwhmz==3,ULM.parameters.NLocalMax = 2;
    else,ULM.parameters.NLocalMax = 3;
    end
end

%% 1 PREPARE INTENSITY MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MATRIX CROPPING %%
% Creates smaller matrix MatInReduced to avoid boundaries, where microbubbles cannot be localized. We avoided padding because padding would
% result in erroneous localization in the boundaries.
MatInReduced = zeros(height,width,numberOfFrames,typename);
MatInReduced(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1,1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:) = MatIn(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1, 1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:);
[height,width,numberOfFrames] = size(MatInReduced);


%% 2 DETECTION AND SELECTION OF MICROBUBBLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DETECTION OF LOCAL MAXIMA %%
% Concatenates MatInReduced into a 2D matrix with time and space in its rows and only space in columns to apply imregionalmax function
% The imregionalmax connectivity is set to default (8), to consider the 8 adjacent pixels in the horizontal(2), vertical(2), and diagonal directions(4).
% Generates an IntensityMatrix (3D) with only local maximal pixels with associated values.
Mat2D = permute(MatInReduced, [1,3,2]); %so that all the frames are in columns
Mat2D = reshape(Mat2D,height*numberOfFrames,width);% Concatenate Matrix
mask2D = imregionalmax(Mat2D); clear Mat2D  % Perform imregionalmax
mask = reshape(mask2D,height,numberOfFrames,width);clear mask2D % reshape concatenated mask
mask = permute(mask,[1,3,2]); % so that we restore (z,x,t) table

IntensityMatrix = MatInReduced.*mask; %Values of intensities at regional maxima

% SELECTION OF MICROBUBBLES %%
% Only the first numberOfParticles highest local max will be kept for localization.
% Other local max will be considered as noise.
% Sort intensites in each frames, and store pixel coordinates
% At the end of this section, spatial and temporal coordinates microbubbles are
% stored into: index_mask_z, index_mask_x, index_numberOfFrames
[tempMatrix,~] = sort(reshape(IntensityMatrix,[],size(IntensityMatrix,3)),1,'descend');

% Remove the last kept intensity values to each frame. This means that you cannot fix an intensity threshold,
% we rely on number of particles. This is key for transparency/parallelization.
IntensityFinal = IntensityMatrix - ones(size(IntensityMatrix)) .* reshape(tempMatrix( ULM.numberOfParticles+1,:),[1 1 numberOfFrames]);
clear tempMatrix
% Construction of the final mask with only the kept microbubbles low resolved and their associated intensity
MaskFinal = (mask.*IntensityFinal)>0;
MaskFinal(isnan(MaskFinal))=0;
MaskFinal = (MaskFinal>0).*IntensityMatrix;

% Preparing intensities and coordinates for further calculation of average, intensities etc...
index_mask  = find(MaskFinal);
[index_mask_z,index_mask_x,index_numberOfFrames]=ind2sub([height, width, numberOfFrames], index_mask);
clear mask IntensityFinal MaskFinal IntensityMatrix
clear index_mask


%% 3 SUBWALENGTH LOCALIZATION OF MICROBUBBLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOCALIZATION OF MICROBUBBLES %%
% The position of each microbubble will be localized with a subwavelength precision.
% For example a microbubble at pixel coordinates [34 67] will be localized at [34.4 66.8]

% Initialize variables averageXc, averageZc which are the super-resolved position of microbubbles
averageXc = nan(1,size(index_mask_z,1),typename);
averageZc = nan(1,size(index_mask_z,1),typename);

for iscat=1:size(index_mask_z,1)
    % For each microbubble, create a 2D intensity matrix of the Region of interest defined by fwhm
    IntensityRoi = MatIn(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat));

    % NLocal max
    % If there are too many localmax in the region of interest, the microbubble shape will be affected and the localization distorted.
    % In that case, we set averageZc, averageXc to NaN value.
    if nnz(imregionalmax(IntensityRoi))>ULM.parameters.NLocalMax
        continue
    end

    % Apply the localization method selected
    % functions are detailed at the end of the code (excepted LocRadialSym which requires an additional function)
    switch lower(ULM.LocMethod)
        case 'radial'
            [Zc,Xc,sigma] = LocRadialSym(IntensityRoi,fwhmz,fwhmx);
        case 'wa'
            [Zc,Xc,sigma] = LocWeightedAverage(IntensityRoi,vectfwhmz,vectfwhmx);
        case 'interp'
            [Zc,Xc,sigma] = LocInterp(IntensityRoi,ULM.parameters.InterpMethod,vectfwhmz,vectfwhmx);
        case 'curvefitting'
            [Zc,Xc,sigma] = curveFitting(IntensityRoi,vectfwhmz,vectfwhmx);
        case 'nolocalization'
            [Zc,Xc,sigma] = NoLocalization(IntensityRoi);
        otherwise
            error('Wrong LocMethod selected')
    end

    % Store the final super-resolved position of the microbubble as its pixel position and an axial/lateral sub-pixel shift.
    averageZc(iscat) = Zc + index_mask_z(iscat);
    averageXc(iscat) = Xc + index_mask_x(iscat);

    % Additional safeguards
    % sigma evaluates the size of the microbubble. If it appears to be too large, the microbubble can be removed (optional)
    if or(sigma<0,sigma>25)
%         averageZc(iscat)=nan;
%         averageXc(iscat)=nan;
%         continue
    end

    % If the final axial/lateral shift is higher that the fwhmz,
    % localization has diverged and the microbubble is ignored.
    if or(abs(Zc)>fwhmz/2,abs(Xc)>fwhmx/2)
        averageZc(iscat)=nan;
        averageXc(iscat)=nan;
        continue
    end
end
keepIndex = ~isnan(averageXc);

ind = sub2ind([height,width,numberOfFrames],index_mask_z(keepIndex),index_mask_x(keepIndex),index_numberOfFrames(keepIndex));
clear index_mask_z index_mask_x IntensityRoi

%% BUILD MATTRACKING %%
% Creating the table which stores the high resolved microbubbles coordinates and the density value
MatTracking = zeros(nnz(keepIndex),4,typename);

MatTracking(:,1) = MatInReduced(ind);       % Initial intensity of the microbubble
MatTracking(:,2) = averageZc(keepIndex);    % Super-resolved axial coordinate
MatTracking(:,3) = averageXc(keepIndex);    % Super-resolved lateral coordinate
MatTracking(:,4) = index_numberOfFrames(keepIndex); % Frame number of the microbubble
clear averageXc averageZc index_numberOfFrames MatInReduced
end

%% ADDITIONAL LOCALIZATION FUNCTIONS

function sigma = ComputeSigmaScat(Iin,Zc,Xc)
%% This function will calculate the Gaussian width of the presupposed peak in the intensity, which we set as an estimate of the width of the microbubble
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
    [Zc,Xc] = localizeRadialSymmetry(Iin,fwhm_z,fwhm_x);
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocRadialSym_sg(Iin,fwhm_z,fwhm_x)
%% function [Zc,Xc,sigma] = LocRadialSym(Iin,fwhm_z,fwhm_x)
    [Zc,Xc] = localizeRadialSymmetry_sg(Iin,fwhm_z,fwhm_x);
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = NoLocalization(Iin)
%% function [Zc,Xc,sigma] = NoLocalization(Iin)
    % position is kept in the center
    Xc = 0;Zc = 0;
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
    % weighted average axial and lateral calculation
    Zc = sum(sum(Iin.*vectfwhm_z',1),2)./sum(Iin(:));
    Xc = sum(sum(Iin.*vectfwhm_x,1),2)./sum(Iin(:));
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)
    % Interpolated base model, the final position if the coordinates of the max intensity interpolated pixel.
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

    Zc = vectfwhm_z(1) -.5 + iz./10 -.05;
    Xc = vectfwhm_x(1) -.5 + ix./10 -.05;
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
% The ROI intensity is fitted with a theorical microbubble model.
    [meshX,meshZ] = meshgrid(vectfwhm_x,vectfwhm_z);
    meshIn = cat(3,meshX,meshZ);

    sigGauss_z = vectfwhm_z(end)*0+1;
    sigGauss_x = vectfwhm_x(end)*0+1;
    myGaussFunc = @(x_pos,mesh_pos)( exp(-(mesh_pos(:,:,1)-x_pos(1)).^2./(2*sigGauss_z^2) - (mesh_pos(:,:,2)-x_pos(2)).^2./(2*sigGauss_x^2)));
    OPTIONS = optimoptions('lsqcurvefit','StepTolerance',.01,'MaxIterations',5,'Display','off');

    % Gaussian fitting
    x_out = lsqcurvefit(myGaussFunc,[0 0],meshIn,double(Iin./max(Iin(:))),[],[],OPTIONS);
    Zc = x_out(2);
    Xc = x_out(1);
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end
