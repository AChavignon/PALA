% radialcenter.m
%
% Created by Baptiste Heiles on 05/09/18
% Inspired from Raghuveer Parthasarathy, The University of Oregon
%%
% Calculates the center of a 2D intensity distribution.
% Method: According to Parthasarathy 2012,uses the radial symmetry of the
% PSF. The gradient of a function that has perfect radial symmetry will
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
%

function [xc,zc,sigma] = localizeRadialSymmetry(I,fwhmz,fwhmx)

% Number of grid points
[Nz,Nx] = size(I);

if nnz(imregionalmax(I))>3 % default 3
    I_interp=interp2(I);
    [Nz_interp,Nx_interp]=size(I_interp);
    fwhmz_interp=fwhmz*2-1;fwhmx_interp=fwhmx*2-1;
    vectfwhmz_interp=-1*round(fwhmz_interp/2):round(fwhmz_interp/2);
    vectfwhmx_interp=-1*round(fwhmx_interp/2):round(fwhmx_interp/2);
    Zc=sum(sum(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp).*vectfwhmz_interp',1),2)./sum(reshape(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp),[],1));
    Xc=sum(sum(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp).*vectfwhmx_interp,1),2)./sum(reshape(I_interp(round(Nz_interp/2)+vectfwhmz_interp,round(Nx_interp/2)+vectfwhmx_interp),[],1));
    zc=Zc/2;xc=Xc/2;
    
else
    % grid coordinates are -n:n, where Nz (or Nx) = 2*n+1
    % grid midpoint coordinates are -n+0.5:n-0.5;
    zm_onerow = (-(Nz-1)/2.0+0.5:(Nz-1)/2.0-0.5)';
    zm = zm_onerow(:,ones(Nx-1, 1));
    xm_onecol = (-(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5);  % Note that y increases "downward"
    xm = xm_onecol(ones(Nz-1, 1),:);
    
    % Calculate derivatives along 45-degree shifted coordinates (u and v)
    % Note that y increases "downward" (increasing row number) -- we'll deal
    % with this when calculating "m" below.
    dIdu = I(1:Nz-1,2:Nx)-I(2:Nz,1:Nx-1);
    dIdv = I(1:Nz-1,1:Nx-1)-I(2:Nz,2:Nx);
    
    % Smoothing --
    h = ones(3)/9;  % simple 3x3 averaging filter
    fdu = conv2(dIdu, h, 'same');
    fdv = conv2(dIdv, h, 'same');
    dImag2 = fdu.*fdu + fdv.*fdv; % gradient magnitude, squared
    
    % Slope of the gradient .  Note that we need a 45 degree rotation of
    % the u,v components to express the slope in the x-y coordinate system.
    % The negative sign "flips" the array to account for y increasing
    % "downward"
    m = -(fdv + fdu) ./ (fdu-fdv);
    
    % *Very* rarely, m might be NaN if (fdv + fdu) and (fdv - fdu) are both
    % zero.  In this case, replace with the un-smoothed gradient.
    NNanm = sum(isnan(m(:)));
    if NNanm > 0
        unsmoothm = (dIdv + dIdu) ./ (dIdu-dIdv);
        m(isnan(m))=unsmoothm(isnan(m));
    end
    % If it's still NaN, replace with zero. (Very unlikely.)
    NNanm = sum(isnan(m(:)));
    if NNanm > 0
        m(isnan(m))=0;
    end
    
    %
    % Almost as rarely, an element of m can be infinite if the smoothed u and v
    % derivatives are identical.  To avoid NaNs later, replace these with some
    % large number -- 10x the largest non-infinite slope.  The sign of the
    % infinity doesn't matter
    try
        m(isinf(m))=10*max(m(~isinf(m)));
    catch
        % if this fails, it's because all the elements are infinite.  Replace
        % with the unsmoothed derivative.  There's probably a more elegant way
        % to do this.
        m = (dIdv + dIdu) ./ (dIdu-dIdv);
    end
    
    
    % Shorthand "b", which also happens to be the
    % y intercept of the line of slope m that goes through each grid midpoint
    b = zm - m.*xm;
    
    % Weighting: weight by square of gradient magnitude and inverse
    % distance to gradient intensity centroid.
    sdI2 = sum(dImag2(:));
    zcentroid = sum(sum(dImag2.*zm))/sdI2;
    xcentroid = sum(sum(dImag2.*xm))/sdI2;
    w  = dImag2./sqrt((zm-zcentroid).*(zm-zcentroid)+(xm-xcentroid).*(xm-xcentroid));
    
    % least-squares minimization to determine the translated coordinate
    % system origin (xc, yc) such that lines y = mx+b have
    % the minimal total distance^2 to the origin:
    % See function lsradialcenterfit (below)
    [zc,xc] = lsradialcenterfit(m, b, w);
    
end

%% Calculate sigma
% Return output relative to upper left coordinate
%xc = xc + (Nx+1)/2.0;%BH
%yc = yc + (Ny+1)/2.0;%BH

% A rough measure of the particle width.
% Not at all connected to center determination, but may be useful for tracking applications;
% could eliminate for (very slightly) greater speed
Isub = I - min(I(:));
[px,pz] = meshgrid(1:Nx,1:Nz);
zoffset = pz - zc+(Nz)/2.0;%BH xoffset = px - xc;
xoffset = px - xc+(Nx)/2.0;%BH yoffset = py - yc;
r2 = zoffset.*zoffset + xoffset.*xoffset;
sigma = sqrt(sum(sum(Isub.*r2))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width

%%

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
% % % to compare between this method and weighted average
% % % if nnz(imregionalmax(I))>1
% % % fwhmz=3;fwhmx=3;
% % % vectfwhmz=-1*round(fwhmz/2):round(fwhmz/2);
% % % vectfwhmx=-1*round(fwhmx/2):round(fwhmx/2);
% % % xt=3 + sum(sum(I.*vectfwhmz',1),2)./sum(I(:));
% % % yt=3 + sum(sum(I.*vectfwhmz',1),2)./sum(I(:));
% % % imagesc(I);hold on; plot( xc+Nx/2, yc+Ny/2,'rx'); plot( xt, yt,'ko');drawnow;hold off;pause
% % % end
end