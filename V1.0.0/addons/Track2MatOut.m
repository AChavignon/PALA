function MatOut = Track2MatOut(Tracks,sizeOut)
%% function MatOut = Track2MatOut(Tracks,sizeOut)
% This function returns a density image from all tracks. Each track is projected on an
% grid image. If the track cross a pixel, the instensity of this pixel is increased by 1.
% At the end, the image represents the number of tracks that has crossed pixels.
% INPUTS:
%       - Tracks is the list of tracks converted in superresolved pixel. The conversion
%       has to be done before sending tracks to this function.
%       - sizeOut is the size of the grid where tracks are projected.
% OUTPUTS:
%       - MatOut is the density map of size sizeOut (generaly 10x10 bigger the size of IQ)
% Created by Arthur Chavignon, 11/04/19
%
% DATE 2020.03.30 - VERSION 1.0.0
% AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
% Code Available under Creative Commons Non-Commercial 4.0
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018 


%% Inputs position are in pixel (for example z=1.2,x=3.54), then converted in index
% in the grid by rounding value (z=1,x=4). The pixel of position (1,4) will be increment
% by 1.

%% select mode
if iscell(Tracks)
    mode = '2D_tracks';
else
    mode = '2D_allin';
end

%% Compute MatOut
size_z = sizeOut(1);
size_x = sizeOut(2);
MatOut=zeros(size_z,size_x); % Initialization of the density image

switch mode
    case '2D_allin'
        %% 2D case - 1 cell
        disp(['Building 2D ULM image (size : ' num2str(sizeOut) ').'])
        % round pixel position into [z,x] index           
        pos_z_round = round(Tracks(:,1));
        pos_x_round = round(Tracks(:,2));

        % remove out of grid bubbles (ie. the grid is too small)
        outP = zeros(numel(pos_z_round),1);
        outP = or(outP,pos_z_round<1);outP = or(outP,pos_z_round>sizeOut(1));
        outP = or(outP,pos_x_round<1);outP = or(outP,pos_x_round>sizeOut(2));

        ind=(sub2ind([size_z size_x],pos_z_round(~outP),pos_x_round(~outP)));

        for j=1:size(ind,1)
            MatOut(ind(j))=MatOut(ind(j))+1;
        end

    case '2D_tracks'
        %% 2D case - 1 cell per tracks
        disp(['Building 2D ULM image (size : ' num2str(sizeOut) ').'])
        for itrack = 1:numel(Tracks)
            % round pixel position into [z,x] index           
            pos_z_round = round(Tracks{itrack}(:,1));  
            pos_x_round = round(Tracks{itrack}(:,2));

            % remove out of grid bubbles (ie. the grid is too small)
            outP = zeros(numel(pos_z_round),1);
            outP = or(outP,pos_z_round<1);outP = or(outP,pos_z_round>sizeOut(1));
            outP = or(outP,pos_x_round<1);outP = or(outP,pos_x_round>sizeOut(2));

            ind=(sub2ind([size_z size_x],pos_z_round(~outP),pos_x_round(~outP)));
            ind = unique(ind);
            
            for j=1:size(ind,1)
                MatOut(ind(j))=MatOut(ind(j))+1;
            end
        end
    
    otherwise
        error('Wrong size selected')
end

end
