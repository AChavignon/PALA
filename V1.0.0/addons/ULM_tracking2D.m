function varargout = ULM_tracking2D(MatTracking,numberOfFrames,FrameRate,max_linking_distance,max_gap_closing,min_length,varargin)
%% function varargout = ULM_tracking2D(MatTracking,numberOfFrames,FrameRate,max_linking_distance,max_gap_closing,min_length,varargin)
% Takes bubbles positions and returns a lists of tracks. The pairing of bubbles is done
% 'simpletracker' (from Jean-Yves Tinevez). It takes a cellstruct conainting the bubbles'
% positions of the current frmae (cell(ii) = positions of bubbles at frame ii).
% Short tracks are remove because not consistent. Long tracks are then interpolated to
% later correctly fill the density Matout. 
% All parameters must have the same unit (bubbles position, max_linking_distance), it can be mm, wavelength, um...
% with 
% Tracking 2D function created by Baptiste Heiles to do 2D tracking
% Different modes have been inplemented: 
%       - nointerp : raw tracks [z x iFrame]
%       - interp : interpolation of track [zi xi] using interp_factor (1/interp_factor more points)
%       - velocityinterp : interpolation of tracks and calculate lateral/axial velocities [z x vz vx timeline]
%       - fightclub : return raw tracks AND interpolated tracks
%
%
% INPUTS:
%       - MatTracking is a[4 Nparticles*Nframes] matrix sorting bubble's information [intensity, z position, x position, frame number]
%       - numberOfFrames
%       - max_linking_distance is the maximal distance between 2 bubbles that can be paired
%       - max_gap_closing is maximal gap of frame accepted to pair bubbles.
%       - min_length is the minimal length of tracks
%       - mode is the method used to postprocess tracks (ie. interpolation, velocity measurement...)
% OUPUTS:
%       - Tracks is a cellstruct containing the list of coordinates of each tracks.
%
% DATE 2020.03.30 - VERSION 1.0.0
% AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
% Code Available under Creative Commons Non-Commercial 4.0
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018 

interp_factor=0.05;
if nargin>6
    mode=lower(varargin{1});
else
    mode='velocityinterp';
end
scale_t=1/FrameRate; %[usec]

%% renormalize to take into account that not all matrices start with frame number 1
minFrame = min(MatTracking(:,4));

MatTracking(:,4) = MatTracking(:,4) - minFrame+1;
index_frames=arrayfun(@(i) find(MatTracking(:,4)==i),[1:numberOfFrames],'UniformOutput',false);

Points=arrayfun(@(i) [MatTracking(index_frames{i},2),MatTracking(index_frames{i},3)],[1:numberOfFrames],'UniformOutput',false);
debug=false;
[ tracks,adjacency_tracks ] = simpletracker(Points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);
n_track=numel(tracks);
all_points = vertcat(Points{:});

count=1;
for i_track = 1 :n_track
    track_id = adjacency_tracks{i_track};
    idFrame = MatTracking(track_id,4);
    track_points = cat(2,all_points(track_id,:),idFrame);
    if length(track_points(:,1))>min_length
        Track{count}=track_points;
        count=count+1;
    end
end

if count==1, disp(['Was not able to find tracks at ',num2str(minFrame)]);Track{1}=[0,0,0,0];end

%% Post porcessing of tracks
Tracks = {};
switch mode            
    case 'nointerp'  
        %% without interpolation, raw tracks
        for i_track = 1:size(Track,2)
            track_points=double(Track{1,i_track});
            
            xi=track_points(:,2);
            zi=track_points(:,1);
            iFrame = track_points(:,3);
                        
            if length(zu)>(min_length)
                Tracks{i_track,1}=cat(2,zi,xi,iFrame);
            end
        end
        
    case 'interp'  
        %% without interpolation, raw tracks
        for i_track = 1:size(Track,2)
            track_points=double(Track{1,i_track});
            
            xi=track_points(:,2);
            zi=track_points(:,1);
            zu = interp1(1:length(zi),smooth(zi,20),1:interp_factor:length(zi));
            xu = interp1(1:length(xi),smooth(xi,20),1:interp_factor:length(xi));
            
            if length(zu)>(min_length)
                Tracks{i_track,1}=cat(2,zu,xu);
            end
        end
        
    case 'velocityinterp'  
        %% without interpolation, raw tracks
        for i_track = 1:size(Track,2)
            track_points=double(Track{1,i_track});
            
            xi=track_points(:,2);
            zi=track_points(:,1);
            TimeAbs= (0:(length(z)-1))*scale_t;

            % Interpolation
            zu = interp1(1:length(zi),smooth(zi,20),1:interp_factor:length(zi));
            xu = interp1(1:length(xi),smooth(xi,20),1:interp_factor:length(xi));
            TimeAbs_interp = interp1(1:length(TimeAbs),TimeAbs,1:interp_factor:length(TimeAbs));
            
            % Velocity
            vzu=diff(zu)./diff(TimeInterp);vzu=[vzu(1);vzu];
            vxu=diff(xu)./diff(TimeInterp);vxu=[vxu(1);vxu];
            
            if length(zu)>(min_length)
                Tracks{i_track,1}=cat(2,zu,xu,vzu,vxu,TimeInterp); %position / velocity / timeline
            end
        end
        
    case 'fightclub'
        %% with and without interpolation
        for i_track = 1:size(Track,2)
            track_points=double(Track{1,i_track});
            
            xi=track_points(:,2);
            zi=track_points(:,1);
            iFrame = track_points(:,3);
            
            if length(zi)>(min_length)
                % store in Tracks position and frame number, used to compare with
                % simulation dataset where absolution positions are available.
                Tracks{i_track,1}=cat(2,zi,xi,iFrame);
            end
            
            % Interpolate tracks for density rendering
            zu = interp1(1:length(zi),smooth(zi,20),1:interp_factor:length(zi));
            xu = interp1(1:length(xi),smooth(xi,20),1:interp_factor:length(xi));
            
%             [zu,izu,~]=unique(round(z*20),'stable');
%             [xu,ixu,~]=unique(round(x*20),'stable');
%             ifin=union(izu,ixu,'stable');
            
%             zu=(z(ifin));xu=(x(ifin));
            if length(zu)>(min_length/interp_factor)
                Tracks_interp{i_track,1}=cat(2,zu',xu',repmat(length(zu).*interp_factor./scale_t.*1e3,size(zu))');
            end
            
        end
        if nargout>1
            varargout{2} = Tracks_interp;
        end
end

Tracks = Tracks(~cellfun('isempty',Tracks));

varargout{1}=Tracks;

end

