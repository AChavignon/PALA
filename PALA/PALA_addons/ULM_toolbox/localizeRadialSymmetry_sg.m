function [zc,xc] = localizeRadialSymmetry_sg(I,fwhmz,fwhmx)
%% function [zc,xc] = localizeRadialSymmetry(I,fwhmz,fwhmx)
% localizeRadiaSymmetry.m
%
% Created by Baptiste Heiles on 05/09/18
% Inspired from Raghuveer Parthasarathy, The University of Oregon
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
% Calculates the center of a 2D intensity distribution.
% Method: The gradient of a function that has perfect radial symmetry will
% point towards the origin. Thus we take the local gradient and construct
% lines through any point with orientation parallel to the local gradient.
% The origin is the point that will minimize the distance between itself
% and all such lines.
%
% Inputs
%   I  : 2D intensity distribution
%        Size need not be an odd number of pixels along each dimension
%   fwhmz, fwhmx : full width at half maximum in direction z and x
% Outputs
%   zc, xc : the center of radial symmetry,
%            px, from center
%   sigma  : Rough measure of the width of the distribution (sqrt. of the
%            second moment of I - min(I));
%            Not determined by the fit -- output mainly for consistency of
%            formatting compared to my other fitting functions, and to get
%            an estimate of the particle "width"

%% Number of grid points
[Nz,Nx] = size(I);

%% Depending on the number of local maxima present in the I window, the algorithm will either do Weighted average localization (if there are too many maxima) or radial symmetry
if nnz(imregionalmax(I))>3 % default 3
%     warning('interp')
    I_interp=interp2(I);% to have a better weighted average localization, we will interpolate (linear interpolation) the I window  on a refined grid formed by dividing the interval between sample values once in each dimension
    [Nz_interp,Nx_interp]=size(I_interp);% Get sizes of the window
    fwhmz_interp=fwhmz*2-1;fwhmx_interp=fwhmx*2-1;% Recalculate the fwhms to match with the interpolation
    vectfwhmz_interp=-1*round(fwhmz_interp/2):round(fwhmz_interp/2);% Calculate the weights in Z direction
    vectfwhmx_interp=-1*round(fwhmx_interp/2):round(fwhmx_interp/2);% Calculate the weights in X direction
    % Calculate the assumed center with the Weighted Average method
    Zc=sum(sum(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp).*vectfwhmz_interp',1),2)./sum(reshape(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp),[],1));
    Xc=sum(sum(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp).*vectfwhmx_interp,1),2)./sum(reshape(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp),[],1));
    zc=Zc/2;xc=Xc/2;% Calculate the actual position on the non-refined grid
    return
end

%% Radial symmetry algorithm

% grid coordinates are -n:n, where Nz (or Nx) = 2*n+1
% grid midpoint coordinates are -n+0.5:n-0.5. Note that z increases "downward"
zm_onerow = (-(Nz-1)/2.0+0.5:(Nz-1)/2.0-0.5)';
zm = zm_onerow(:,ones(Nx-1, 1));
xm_onecol = (-(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5);
xm = xm_onecol(ones(Nz-1, 1),:);

% Calculate derivatives along 45-degree shifted coordinates (u and v) Please refer to Appendix 2 of the publication attached to this code for basis definition
dIdu = I(1:Nz-1,2:Nx)-I(2:Nz,1:Nx-1);% Gradient along the u vector
dIdv = I(1:Nz-1,1:Nx-1)-I(2:Nz,2:Nx);% Gradient along the v vector

% Smoothing the gradient of the I window
h = ones(3)/9;
fdu = conv2(dIdu, h, 'same');% Convolution of the gradient with a simple averaging filter
fdv = conv2(dIdv, h, 'same');
dImag2 = fdu.*fdu + fdv.*fdv; % Squared gradient magnitude

% Slope of the gradient . Please refer to appendix 2 of the publication attached to this code for basis/orientation
m = -(fdv + fdu) ./ (fdu-fdv);

% Check if m is NaN (which can happen when fdu=fdv). In this case, replace with the un-smoothed gradient.
NNanm = sum(isnan(m(:)));
if NNanm > 0
    unsmoothm = (dIdv + dIdu) ./ (dIdu-dIdv);
    m(isnan(m))=unsmoothm(isnan(m));
end

% If it's still NaN, replace with zero and we'll deal with this later
NNanm = sum(isnan(m(:)));
if NNanm > 0
    m(isnan(m))=0;
end

% Check if m is inf (which can happen when fdu=fdv).
try
    m(isinf(m))=10*max(m(~isinf(m)));
catch
    % Replace m with the unsmoothed gradient
    m = (dIdv + dIdu) ./ (dIdu-dIdv);
end

% Calculate the z intercept of the line of slope m that goes through each grid midpoint
b = zm - m.*xm;

% Weight the intensity by square of gradient magnitude and inverse
% distance to gradient intensity centroid. This will increase the intensity of areas close to the initial guess
sdI2 = sum(dImag2(:));
zcentroid = sum(sum(dImag2.*zm))/sdI2;% Initial guess of the centroid in z
xcentroid = sum(sum(dImag2.*xm))/sdI2;% Initial guess of the centroid in x
w  = dImag2./sqrt((zm-zcentroid).*(zm-zcentroid)+(xm-xcentroid).*(xm-xcentroid));

% least-squares minimization to determine the translated coordinate
% system origin (xc, yc) such that lines y = mx+b have
% the minimal total distance^2 to the origin:
% See function lsradialcenterfit (below)
[zc,xc] = lsradialcenterfit(m, b, w);

end

% We'll code the least square solution function separately as we could find the solution with another implementation
function [zc,xc] = lsradialcenterfit(m, b, w)
    % least squares solution to determine the radial symmetry center

    % inputs m, b, w are defined on a grid
    % w are the weights for each point
    wm2p1 = w./(m.*m+1);
    sw  = sum(sum(wm2p1));
    smmw = sum(sum(m.*m.*wm2p1));
    smw  = sum(sum(m.*wm2p1));
    smbw = sum(sum(m.*b.*wm2p1));
    sbw  = sum(sum(b.*wm2p1));
    det = smw*smw - smmw*sw;
    xc = (smbw*sw - smw*sbw)/det;    % relative to image center
    zc = (smbw*smw - smmw*sbw)/det; % relative to image center
end
