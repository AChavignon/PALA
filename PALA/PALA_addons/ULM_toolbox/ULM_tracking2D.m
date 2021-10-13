function varargout = ULM_tracking2D(MatTracking,ULM,varargin)
%% function varargout = ULM_tracking2D(MatTracking,ULM,varargin)
% Takes bubbles positions and returns a lists of tracks. The pairing of bubbles is based
% on partial assignment and minimization of the total distance with 'simpletracker' and
% munkres implementation (from Jean-Yves Tinevez https://github.com/tinevez/simpletracker).
% It takes a cellstruct containing the bubbles' positions of the current frame (cell(ii) = positions of bubbles at frame ii).
% Short tracks are remove because not consistent. Long tracks are then interpolated to
% later correctly fill the density MatOut.
% All parameters must have the same unit (bubbles position, max_linking_distance), it can be mm, wavelength, um...
%
% Different modes have been implemented:
%       - nointerp : raw tracks [z x iFrame]
%       - interp : interpolation of track [zi xi] using interp_factor (1/interp_factor more points)
%       - velocityinterp : interpolation of tracks and calculate lateral/axial velocities [z x vz vx timeline]
%       - pala : return raw tracks AND interpolated tracks
%
% INPUTS:
%       - MatTracking is a[4 Nparticles*Nframes] matrix sorting bubble's information [intensity, z position, x position, frame number]
%       - numberOfFrames
%       - ULM.max_linking_distance is the maximal distance between 2 bubbles that can be paired
%       - ULM.max_gap_closing is maximal gap of frame accepted to pair bubbles.
%       - ULM.min_length is the minimal length of tracks
%       - ULM.scale is the scale of data, [scale_z scale_x scale_t] (scale_z=scale_x=1, scale_t=1/FrameRate for velocity)
%       - mode is the method used to post process tracks (ie. interpolation, velocity measurement...)
% OUTPUTS:
%       - Tracks is a cellstruct containing the list of coordinates of each tracks.
%
% This function was created by Baptiste Heiles 07/04/17, last modifications Arthur Chavignon, 18/11/19
%
% Requirements: this function requires a two additional functions from the MathWorks
% Community and should be downloaded and add to the path.
%       - 'simpletracker.m' : Jean-Yves Tinevez (2020). simpletracker (https://www.github.com/tinevez/simpletracker), GitHub. Retrieved April 22, 2020.
%           https://fr.mathworks.com/matlabcentral/fileexchange/34040-simpletracker)
%       - 'munkres.m': by Yi Cao at Cranfield University on 17th June 2008 (http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html)
%           https://fr.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

if nargin>2
    mode=lower(varargin{1});
else
    mode='velocityinterp'; % default mode
end
%% Define local parameters
% The interp_factor is required for tracks interpolation. As we want to reconstruct
% rendering with a smaller pixel grid, tracks are interpolated to get at least 1 point in
% the grid. It can be calculated with the maximal distance between two points of a same
% track divided by the resolution factor for rendering (res)
interp_factor = 1/ULM.max_linking_distance/ULM.res*.8;

% We observed that smoothing tracks provide a better rendering and velocity estimation. It
% may correct pairings mistakes, and localization errors along a track. It can be adjust
% by user.
smooth_factor = 20; %
numberOfFrames = ULM.size(3);

%% Renormalizes to take into account that not all matrices start with frame number 1
minFrame = min(MatTracking(:,4));

MatTracking(:,4) = MatTracking(:,4) - minFrame+1;
index_frames=arrayfun(@(i) find(MatTracking(:,4)==i),[1:numberOfFrames],'UniformOutput',false);

Points=arrayfun(@(i) [MatTracking(index_frames{i},2),MatTracking(index_frames{i},3)],[1:numberOfFrames],'UniformOutput',false);
debug=false;
[ Simple_Tracks,Adjacency_Tracks ] = simpletracker(Points,...
    'MaxLinkingDistance', ULM.max_linking_distance, ...
    'MaxGapClosing', ULM.max_gap_closing, ...
    'Debug', debug);
n_tracks=numel(Simple_Tracks);
all_points = vertcat(Points{:});

count=1;Tracks_raw = {};
for i_track = 1:n_tracks
    track_id = Adjacency_Tracks{i_track};
    idFrame = MatTracking(track_id,4);
    track_points = cat(2,all_points(track_id,:),idFrame);
    if length(track_points(:,1))>ULM.min_length
        Tracks_raw{count}=track_points;
        count=count+1;
    end
end

if count==1
    disp(['Was not able to find tracks at ',num2str(minFrame)]);
    Tracks_out{1}=[0,0,0,0];varargout{1}=Tracks_out;
    if nargout>1,varargout{2} = Tracks_out;end
    return
end

%% Post processing of tracks
Tracks_out = {};
switch lower(mode)
    case 'nointerp'
        % without interpolation, raw tracks
        for i_track = 1:size(Tracks_raw,2)
            track_points=double(Tracks_raw{1,i_track});
            xi=track_points(:,2);
            zi=track_points(:,1);
            iFrame=track_points(:,3);

            if length(zi)>ULM.min_length
                Tracks_out{i_track,1}=cat(2,zi,xi,iFrame);
            end
        end

    case 'interp'
        % with tracks interpolation
        for i_track = 1:size(Tracks_raw,2)
            track_points=double(Tracks_raw{1,i_track});
            xi=track_points(:,2);
            zi=track_points(:,1);
            zu=interp1(1:length(zi),smooth(zi,smooth_factor),1:interp_factor:length(zi));
            xu=interp1(1:length(xi),smooth(xi,smooth_factor),1:interp_factor:length(xi));

            if length(zi)>ULM.min_length
                Tracks_out{i_track,1}=cat(2,zu,xu);
            end
        end

    case 'velocityinterp'
        % with tracks interpolation and velocity computation
        for i_track = 1:size(Tracks_raw,2)
            track_points=double(Tracks_raw{1,i_track});
            xi=track_points(:,2);
            zi=track_points(:,1);
            TimeAbs=(0:(length(zi)-1))*ULM.scale(3);

            % Interpolation of spatial and time components
            zu=interp1(1:length(zi),smooth(zi,smooth_factor),1:interp_factor:length(zi));
            xu=interp1(1:length(xi),smooth(xi,smooth_factor),1:interp_factor:length(xi));
            TimeAbs_interp = interp1(1:length(TimeAbs),TimeAbs,1:interp_factor:length(TimeAbs));

            % Velocity
            vzu=diff(zu)./diff(TimeAbs_interp);vzu=[vzu(1),vzu];
            vxu=diff(xu)./diff(TimeAbs_interp);vxu=[vxu(1),vxu];

            if length(zi)>ULM.min_length
                Tracks_out{i_track,1}=single(cat(2,zu',xu',vzu',vxu',TimeAbs_interp')); %position / velocity / timeline
            end
        end

    case {'pala'}
        % with and without interpolation, dedicated to PALA comparison of localization algorithms.
        for i_track = 1:size(Tracks_raw,2)
            track_points=double(Tracks_raw{1,i_track});
            xi=track_points(:,2);
            zi=track_points(:,1);
            iFrame = track_points(:,3);

            if length(zi)>ULM.min_length
                % store in Tracks position and frame number, used to compare with
                % simulation dataset where absolution positions are available.
                Tracks_out{i_track,1}=cat(2,zi,xi,iFrame);
            end

            % Interpolate tracks for density rendering
            zu = interp1(1:length(zi),smooth(zi,smooth_factor),1:interp_factor:length(zi));
            xu = interp1(1:length(xi),smooth(xi,smooth_factor),1:interp_factor:length(xi));
            dd = sqrt(diff(xu).^2+diff(zu).^2); % curvilinear abscissa
            vmean = sum(dd)./(length(zi))/ULM.scale(3); % averaged velocity of the track in [unit]/s

            if length(zi)>ULM.min_length
%                 Tracks_interp{i_track,1}=cat(2,zu',xu',repmat(length(zu).*interp_factor./ULM.scale(3).*1e3,size(zu))');
                Tracks_interp{i_track,1}=cat(2,zu',xu',vmean.*ones(size(zu')));
            end

        end
        if nargout>1
            varargout{2} = Tracks_interp;
        end
end
Tracks_out = Tracks_out(~cellfun('isempty',Tracks_out));
varargout{1}=Tracks_out;
end
